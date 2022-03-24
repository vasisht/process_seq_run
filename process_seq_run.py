#!/usr/bin/python
import sys
from collections import defaultdict
import json
import shlex
import os
import shutil
import subprocess
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import traceback
import time
import glob
from utils import *
import argparse
import ConfigParser
import re


def demultiplex_run(run_path, output_path, settings, umi=False,
                    umi_single_end=False):
    if os.path.isdir(output_path):
        # Make sure we don't delete data above current dir
        # We are assuming the parent directories for run_path and output_path
        # are the same
        parent_dir, run = os.path.split(run_path)
        if output_path.startswith(os.path.abspath(parent_dir) + '/'):
            shutil.rmtree(output_path)
    try:
        umi_opts = ""
        if umi:
            base_mask = 'y*,i8y*,y*'
            if umi_single_end:
                base_mask = 'y*,i8y*'
            umi_opts = (" --use-bases-mask %s "
                        " --minimum-trimmed-read-length=0 "
                        " --mask-short-adapter-reads=0 " % base_mask)
        cmd_opts = "-r 8 -d 8 -p 8 --ignore-missing-bcls --no-lane-splitting"
        cmd = "%s -R %s -o %s %s %s" % (settings['bcl2fastq'], run_path,
                                        output_path, cmd_opts, umi_opts)
        sys.stderr.write("Demultiplexing command\n%s\n" % cmd)
        subprocess.check_call(shlex.split(cmd))
    except Exception as e:
        print traceback.format_exc()
        send_email(settings['from'], settings['to'],
                   settings['server'], settings['password'],
                   settings['port'],
                   run, traceback.format_exc())

    # Generate SAV summary
    summary_file = os.path.join(run_path, "SAV_summary.tsv")
    summary_cmd = "%s %s" % (settings['summary'], run_path)
    sys.stderr.write("Generating summary stats\n%s\n" % summary_cmd)
    with open(summary_file, 'w') as f:
        subprocess.check_call(shlex.split(summary_cmd), stdout=f)

    # Generate index summary
    index_interop = os.path.join(run_path, "InterOp/IndexMetricsOut.bin")
    if os.path.exists(index_interop):
        index_file = os.path.join(run_path, "index_summary.csv")
        index_cmd = "%s %s" % (settings['interopcsv'], index_interop)
        sys.stderr.write("Generating index summary\n%s\n" % index_cmd)
        with open(index_file, 'w') as f:
            subprocess.check_call(shlex.split(index_cmd), stdout=f)
    else:
        sys.stderr.write("Index metrics are not present\n")


def get_dirs_to_process(datadir):
    run_dirs = []
    for dirs in os.listdir(datadir):
        run_dir = os.path.join(datadir, dirs)
        if os.path.isdir(run_dir):
            if os.path.exists(os.path.join(run_dir, 'RTAComplete.txt')):
                if not os.path.exists(os.path.join(run_dir,
                                                   'DemuxComplete.txt')):
                    run_dirs.append(run_dir)
    return run_dirs


def mark_demux_complete(run_dir):
    out = open(run_dir + '/DemuxComplete.txt', 'w')
    out.write('Demultiplexing is complete\n')
    out.close()


def send_email(sender, rcpt, smtp_server, password, port, subject, body):
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = sender
    msg["To"] = rcpt
    msg.attach(MIMEText(body))
    server = smtplib.SMTP(smtp_server, int(port))
    server.starttls()
    server.login(sender, password)
    server.sendmail(msg["From"], msg["To"].split(","), msg.as_string())
    server.quit()


def pretty_print_run_stats(run, overall_metrics, read_summary, lane_summary,
                           index_metrics, output_dir, web_loc):
    run_name = os.path.split(run)[1]
    clustering = "Clustering Stats: \n\n" + \
        "Cluster Density: %.2f k/mm2\n" % overall_metrics['Cluster Density'] + \
        "Cluster Density (PF): %.2f k/mm2\n" % overall_metrics['Cluster Density (PF)'] + \
        "Number of Clusters: %d\n" % overall_metrics['Number of Clusters'] + \
        "Number of Clusters (PF): %d\n" % overall_metrics['Number of Clusters (PF)'] + \
        "Percentage Clusters (PF): %.2f" % overall_metrics['Percentage Clusters (PF)'] + "%\n" + \
        "Percentage PhiX Aligned: %.2f" % overall_metrics['Percentage PhiX Aligned'] + "%\n" + \
        "Error Rate: %.2f" % overall_metrics['Error%'] + "%\n"

    read_fields = ['Level', 'Yield', 'Aligned%', 'Error Rate%', '%>=Q30']
    lane_fields = ['Lane', 'Read', 'Density(k/mm2)', 'Cluster PF%',
                   '%>=Q30', 'Aligned%', 'Error%']
    read_data = "\n\nRead Level Summary: \n\n"
    read_data += print_summary(read_summary, read_fields, 15)

    lane_data = "\n\n Lane Summary \n\n"
    lane_data += print_summary(lane_summary, lane_fields, 17)

    index_data = "\n\n Index Summary \n\n"
    index_fields = ['Sample', 'SampleName', 'Index', 'Counts']
    index_data += print_summary(index_metrics, index_fields, 25)

    out_txt = "\n\n\nOutput directory is %s" % output_dir
    out_txt += "\nhttp://%s/%s/\n" % (web_loc, output_dir)
    body = clustering + read_data + lane_data + index_data + out_txt
    subject = "Processed %s" % run_name
    return subject, body


def process_run(run, settings, upload, nomail, upload_only):
    sender = settings['from']
    rcpt = settings['to']
    smtp_server = settings['server']
    smtp_password = settings['password']
    smtp_port = settings['port']
    samplesheet = os.path.join(run, "SampleSheet.csv")
    sav_summary = os.path.join(run, "SAV_summary.tsv")
    index_summary = os.path.join(run, "index_summary.csv")
    if not os.path.isfile(samplesheet):
        # Check if a csv file is present in the directory
        csvfiles = glob.glob(run + "/*.csv")
        num_csvfiles = len(csvfiles)
        if num_csvfiles == 1:
            shutil.copy(csvfiles[0], samplesheet)
        elif num_csvfiles > 1:
            body = "SampleSheet.csv is absent and too many csv files present"
            send_email(sender, rcpt, smtp_server, smtp_password, smtp_port,
                       run, body)
        elif num_csvfiles == 0:
            body = "No csv files present, cannot demultiplex"
            send_email(sender, rcpt, smtp_server, smtp_password, smtp_port,
                       run, body)
    print "Parsing SampleSheet"
    exp_details = parse_samplesheet(samplesheet)
    if re.search("\%", exp_details["experiment"]):
        # Munich style samplesheet
        prefix, suffix = exp_details["experiment"].split("%")
        dir_path, run_name = os.path.split(run)
        output_prefix = "%s/%s_%s" % (dir_path, prefix, run_name)
        output_suffix = suffix
    else:
        output_prefix = "_".join(run.split("_")[0:-1])
        if re.search("_", exp_details["experiment"]):
            output_suffix = "_".join(exp_details["experiment"].split("_")[1:])
        else:
            output_suffix = exp_details["experiment"]
    output_dir = output_prefix + "_" + output_suffix
    if not upload_only:
        print "Demultiplexing %s" % run
        umi = False
        umi_single_end = False
        if re.search("umi", exp_details["description"].lower()) or re.search("umi", exp_details["experiment"].lower()):
            umi = True
            if len(exp_details['read_lengths']) == 1:
                umi_single_end = True
        demultiplex_run(run, output_dir, settings, umi, umi_single_end)
        if umi:
            # Add barcodes to read names
            processUMI(output_dir, exp_details, umi_single_end)
    print "Parsing SAV Summary"
    read_summary, lane_summary, overall_metrics, index_metrics = summarize_SAV(
        run)
    print "Adding sample names to index"
    index_metrics = add_samplenames_to_index(index_metrics, exp_details)
    exp_details = update_expdetails_with_counts(exp_details, index_metrics)
    exp_details["lane_summary"] = lane_summary
    exp_details["read_summary"] = read_summary
    exp_details["run"] = os.path.split(output_dir)[1]
    exp_details.update(overall_metrics)
    run_json = output_dir + "/run_details.json"
    with open(run_json, "w") as f:
        f.write(json.dumps(exp_details, indent=4, sort_keys=True))
    shutil.copy(samplesheet, output_dir)
    shutil.copy(sav_summary, output_dir)
    shutil.copy(index_summary, output_dir)
    mark_demux_complete(run)
    subject, body = pretty_print_run_stats(run, overall_metrics, read_summary,
                                           lane_summary, index_metrics,
                                           output_dir, settings['web_loc'])
    if not nomail:
        send_email(sender, rcpt, smtp_server, smtp_password, smtp_port,
                   subject, body)
    if upload:
        s4opts = False
        if settings['s4cmd'].lower() == 'true':
            s4opts = True
        error = upload_run_to_S3(settings['s3cfg'], output_dir, settings[
                                 's3folder'], settings['region'], s4opts)
        if error != 0:
            subject = "Upload failure for %s" % output_dir
            body = "Failed to upload %s after 5 tries" % output_dir
            if nomail:
                print body
            else:
                send_email(sender, rcpt, smtp_server, smtp_password, smtp_port,
                           subject, body)
        else:
            # Upload run information to database
            data, exists = create_run_in_db(settings['dbserver'],
                                            settings['dbuser'],
                                            settings['dbpasswd'], run_json)
            if exists:
                subject = 'Run %s already in db' % exp_details['run']
            else:
                subject = "Created %s in db" % exp_details['run']
            send_email(sender, rcpt,  smtp_server, smtp_password,
                       smtp_port, subject, json.dumps(data, indent=4))


def main(settings, upload, nomail, upload_only):
    while 1 == 1:
        run_list = get_dirs_to_process(settings['run_directory'])
        if run_list:
            print "Processing Runs"
            for run in run_list:
                print run
            for run in run_list:
                process_run(run, settings, upload, nomail, upload_only)

        print "Waiting to process runs, sleeping 10 min"
        time.sleep(600)


parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("-c", "--config", help="config file", default='config.cfg')
parser.add_argument("-u", "--upload", help="Upload the data to S3",
                    action="store_true", default=False)
parser.add_argument("--upload_only", help="Don't demux, upload existing data",
                    action="store_true", default=False)
parser.add_argument("--nomail", action="store_true", default=False)
group.add_argument("-i", "--input_dir", help="Run folder to analyze")
group.add_argument(
    "-d", "--daemon", help="Run program as a daemon", action="store_true")

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(0)

args = parser.parse_args()

config = ConfigParser.ConfigParser()
config.read(args.config)
# Set up constants from config
settings = dict(config.items('Data'))
sender = settings['from']
rcpt = settings['to']
smtp_pass = settings['password']
smtp_server = settings['server']
smtp_port = settings['port']


if args.daemon:
    try:
        main(settings, args.upload, args.nomail, args.upload_only)
    except Exception as e:
        print traceback.format_exc()
        if not args.nomail:
            send_email(sender, rcpt, smtp_server, smtp_pass,
                       smtp_port, "Processing error",
                       traceback.format_exc())
elif args.input_dir:
    try:
        process_run(args.input_dir, settings, args.upload, args.nomail,
                    args.upload_only)
    except Exception as e:
        print traceback.format_exc()
        if not args.nomail:
            send_email(sender, rcpt, smtp_server, smtp_pass,
                       smtp_port,  "Processing error",
                       traceback.format_exc())
