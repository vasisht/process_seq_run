#!/usr/bin/python
import re
import datetime
import os
from collections import defaultdict
import sys
import json
import glob
import shutil
from subprocess import check_call
import rnaseq_rest.helpers as rest
from AddUmiNugen import *
from multiprocessing import Pool


def find_string_index(input_str, input_list):
    regex = re.compile(".*" + input_str + ".*")
    item = [m.group(0) for i in input_list for m in [regex.search(i)] if m]
    return input_list.index(item[0])


def parse_samplesheet(infile):

    exp_details = {}
    samples = []
    with open(infile) as f:
        samplesheet = f.read().splitlines()

    header = find_string_index('Header', samplesheet)
    reads = find_string_index('Reads', samplesheet)
    settings = find_string_index('Settings', samplesheet)
    data = find_string_index('Data', samplesheet)

    for lines in samplesheet[header:reads]:
        if re.match("Investigator", lines):
            exp_details["investigator"] = lines.rstrip().split(",")[
                1].replace(" ", "_")
        elif re.match("Experiment", lines):
            exp_details["experiment"] = lines.rstrip().split(",")[
                1].replace(" ", "_")
        elif re.match("Description", lines):
            exp_details["description"] = lines.rstrip().split(",")[
                1].replace(" ", "_")
            if re.search("\%", exp_details["experiment"]):
                # Munich style samplesheet
                exp_details['date'] = datetime.datetime.fromtimestamp(
                    os.path.getctime(infile)).isoformat()
            else:
                try:
                    exp_details['date'] = datetime.datetime.strptime(
                        exp_details["experiment"].split("_")[0],
                        '%Y%b%d').isoformat()
                except ValueError:
                    exp_details['date'] = datetime.datetime.fromtimestamp(
                        os.path.getctime(infile)).isoformat()
        elif re.match("Annotations", lines):
            exp_details["annotations"] = lines.rstrip().split(
                ",")[1].split("|")

    sample_fields = samplesheet[data + 1].split(",")
    read_data = ' '.join(samplesheet[reads + 1:settings])
    read_lengths = re.findall('[0-9]+', read_data)
    exp_details['read_lengths'] = read_lengths

    for lines in samplesheet[data + 2:]:
        if "annotations" in exp_details:
            lines_data = lines.split(",")
            sample_data = dict(zip(sample_fields[:-1], lines_data[:-1]))
            sample_data["annotations"] = dict(
                zip(exp_details["annotations"], lines_data[-1].split("|")))
            samples.append(sample_data)
        else:
            samples.append(dict(zip(sample_fields, lines.split(","))))

    exp_details["samples"] = samples
    return exp_details


def print_summary(summary, fields=[], width=15):
    out_str = ""
    if not fields:
        fields = summary[0].keys()
    for item in fields:
        out_str += "%-*s" % (width, item)
    out_str += "\n"
    for item in summary:
        for val in fields:
            out_str += "%-*s" % (width, item[val])
        out_str += "\n"
    return out_str


def summarize(summary, field, avg=True):
    vals = {}
    for item in summary:
        vals[item['Lane']] = float(item[field].split('+/-')[0])
    if avg:
        return sum(vals.values()) / len(vals.values())
    else:
        return sum(vals.values())


def count_indices(rundir):
    index_counts = defaultdict(int)
    samples = {}
    empty_index = [{'Sample': '', 'Index': '', 'Counts': ''}]
    index_metrics = []
    index_file = os.path.join(rundir, "index_summary.csv")
    if os.path.isfile(index_file):
        infile = open(index_file, 'r')
        for lines in infile:
            if lines[0] == '#' or lines[0] == 'L':
                continue
            else:
                fields = lines.rstrip().split(",")
                samples[fields[4]] = fields[3]
                # Only count for read 1 to prevent double counting
                if fields[2] == '1':
                    index_counts[fields[3]] += int(fields[6])
        infile.close()
        if samples:
            for sample in samples.keys():
                index_metrics.append({'Sample': sample, 'Index': samples[
                                     sample],
                    'Counts': index_counts[samples[sample]]})
            try:
                index_metrics = sorted(
                    index_metrics, key=lambda x: int(x['Sample']))
            except ValueError:
                index_metrics = sorted(
                    index_metrics, key=lambda x: x['Sample'])
            return index_metrics
        else:
            return empty_index
    else:
        return empty_index


def add_samplenames_to_index(index_metrics, exp_details):
    if index_metrics[0]['Sample'] == '':
        index_metrics[0]['SampleName'] = ''
        return index_metrics
    else:
        for index in index_metrics:
            sample = filter(lambda x: x['Sample_ID'] == index[
                            'Sample'], exp_details['samples'])
            if sample:
                index['SampleName'] = sample[0]['Sample_Name']
            else:
                index['SampleName'] = 'Unknown'
        return index_metrics


def summarize_SAV(rundir):
    summary_file = os.path.join(rundir, "SAV_summary.tsv")
    with open(summary_file, 'r') as f:
        data = f.read().splitlines()
    version = find_string_index('# Version', data)
    lane = find_string_index('#Per', data)

    fields = data[version + 1].rstrip().split('\t')
    read_summary = []
    for lines in data[version + 2:lane]:
        read_summary.append(dict(zip(fields, lines.rstrip().split('\t'))))

    fields = data[lane + 1].rstrip().split('\t')
    lane_summary = []
    for lines in data[lane + 2:]:
        lane_summary.append(dict(zip(fields, lines.rstrip().split('\t'))))

    read_fields = ['Level', 'Yield', 'Aligned%', 'Error Rate%', '%>=Q30']
    lane_fields = ['Lane', 'Read', 'Density(k/mm2)', 'Cluster PF%',
                   '%>=Q30', 'Aligned%', 'Error%']

    non_index_reads = filter(
        (lambda x: x['Level'] == 'Non-indexed'), read_summary)[0]
    overall_metrics = {}
    overall_metrics['Cluster Density'] = summarize(
        lane_summary, 'Density(k/mm2)')
    overall_metrics['Cluster Density (PF)'] = summarize(
        lane_summary, 'Density PF')
    overall_metrics['Number of Clusters'] = int(
        summarize(lane_summary, 'Reads', False))
    overall_metrics['Number of Clusters (PF)'] = int(
        summarize(lane_summary, 'Reads PF', False))
    overall_metrics['Percentage Clusters (PF)'] = float(overall_metrics[
        'Number of Clusters (PF)']) * 100.0 / float(overall_metrics[
                                                    'Number of Clusters'])
    overall_metrics['Percentage PhiX Aligned'] = float(
        non_index_reads['Aligned%'])
    overall_metrics['Error%'] = float(non_index_reads['Error Rate%'])

    index_metrics = count_indices(rundir)
    if not index_metrics[0]['Counts'] == '':
        total_indexed_reads = 0
        for indexes in index_metrics:
            total_indexed_reads += indexes['Counts']
        unindexed_reads = overall_metrics[
            'Number of Clusters (PF)'] - total_indexed_reads
        index_metrics.append(
            {'Sample': 'Unknown', 'Index': 'Unknown',
             'Counts': unindexed_reads})

    return read_summary, lane_summary, overall_metrics, index_metrics


def upload_run_to_S3(s3cfg, output_dir, s3_folder, region='us-east-1',
                     s4=False):
    # For Frankfurt use eu-central-1 region
    # TO DO Add ability to specify region to s4cmd
    if s4:
        s4cmd = "s4cmd -p %s dsync --retry 10 -c 12" % s3cfg
        s4cmdopts = (" --ignore-empty-source --retry 5 -c 12 -f -r -s "
                     "--API-ServerSideEncryption=AES256 ")
        run_name = os.path.split(output_dir)[1]
        error = os.system("%s %s %s %s%s" % (
            s4cmd, s4cmdopts, output_dir, s3_folder, run_name))
        # TO DO s4cmd barfs on uploading empty files. Need to fix
        error = 0
    else:
        error = 1
        count = 0
        s3cmd = "s3cmd -c %s sync --limit-rate=10m" % s3cfg
        while (error and count <= 5):
            error = os.system("%s --region=%s --server-side-encryption %s %s" %
                              (s3cmd, region, output_dir, s3_folder))
            count += 1
    return error


def update_expdetails_with_counts(exp_details, index_metrics):
    samples = exp_details['samples']
    for sample in samples:
        sample_data = filter(lambda x: x['SampleName'] == sample[
                             'Sample_Name'], index_metrics)
        if sample_data:
            counts = sample_data[0]['Counts']
        else:
            # Handle cases where the initial samplesheet in incorrect
            # leading to incorrect counts in IndexMetrics.bin
            counts = 0
        sample['Counts'] = counts
    exp_details['samples'] = samples
    return exp_details


def create_run_in_db(dbserver, dbuser, dbpasswd, runjson):
    rest.init(dbserver, username=dbuser, password=dbpasswd)
    with open(runjson, 'r') as f:
        run_data = json.load(f)
    run_exists = rest.queryRNASeqRuns({'name': run_data['run']})
    if run_exists:
        return run_exists, True
    else:
        data = rest.sendRNASeqRunDetail(runjson)
        return data, False

def makedir(path):
    try:
        os.mkdir(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def add_UMI_to_read(indir):
    sys.stderr.write('Processing %s\n' % indir)
    raw_fq = glob.glob(indir + '/*.gz')
    add_umi(indir, indir)
    makedir(indir + '/raw_data')
    for fq in raw_fq:
        sys.stderr.write('moving %s to %s\n' % (fq, indir + '/rawdata'))
        shutil.move(fq, indir + '/raw_data')


def add_UMI_to_read_se(indir):
    sys.stderr.write('Processing %s\n' % indir)
    raw_fq = glob.glob(indir + '/*.gz')
    add_umi_se(indir, indir)
    makedir(indir + '/raw_data')
    for fq in raw_fq:
        sys.stderr.write('moving %s to %s\n' % (fq, indir + '/rawdata'))
        shutil.move(fq, indir + '/raw_data')


def processUMI(rundir, run_details, single_end=False):
    project_name = run_details['samples'][0]['Sample_Project']
    sample_list = os.listdir(rundir + '/' + project_name)
    project_dir = rundir + '/' + project_name
    sys.stderr.write('Base directory is %s\n' % project_dir)
    sample_dirs = [project_dir + '/' + i for i in sample_list]
    pool = Pool(processes=8)
    if single_end:
        pool.map(add_UMI_to_read_se, sample_dirs)
    else:
        pool.map(add_UMI_to_read, sample_dirs)
