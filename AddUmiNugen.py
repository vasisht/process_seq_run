#!/usr/bin/python
import gzip
import glob
import itertools
import io
import subprocess
from readfq import readfq


def writefq(fh, data):
    fh.write('@' + data[0] + '\n')
    fh.write(data[1] + '\n')
    fh.write('+\n')
    fh.write(data[2] + '\n')


def add_umi(input_dir, output_dir):
    R1 = glob.glob(input_dir + '/*R1*fastq*')[0]
    R2 = glob.glob(input_dir + '/*R3*fastq*')[0]
    umi = glob.glob(input_dir + '/*R2*fastq*')[0]
    sample_file = R1.split('/')[-1]
    sample_file_split =  sample_file.split('_')
    if len(sample_file_split) == 4:
        sample_name = sample_file_split[0]
    else:
        sample_name = '_'.join(sample_file_split[:-3])
    r1 = open(output_dir + '/' + sample_name + '_R1_UMI_001.fastq.gz', 'w')
    r2 = open(output_dir + '/' + sample_name + '_R2_UMI_001.fastq.gz', 'w')
    # Use subprocess to pipe to gzip
    r1_out = subprocess.Popen('gzip', stdin=subprocess.PIPE, stdout=r1)
    r2_out = subprocess.Popen('gzip', stdin=subprocess.PIPE, stdout=r2)

    with io.BufferedReader(gzip.open(R1, 'rb')) as read1, \
            io.BufferedReader(gzip.open(R2, 'rb')) as read2, \
            io.BufferedReader(gzip.open(umi, 'rb')) as uid:
        for data1, data2, umidata in itertools.izip(readfq(read1),
                                                    readfq(read2),
                                                    readfq(uid)):
            r1_name, r1_suffix = data1[0].split()
            r2_name, r2_suffix = data2[0].split()
            umi_name = umidata[0].split()[0]
            umi_id = umidata[1]
            assert(r1_name == r2_name == umi_name), 'Mismatch in read names'
            r1_umi = r1_name + ':' + umi_id + ' ' + r1_suffix
            r2_umi = r2_name + ':' + umi_id + ' ' + r2_suffix
            writefq(r1_out.stdin, (r1_umi, data1[1], data1[2]))
            writefq(r2_out.stdin, (r2_umi, data2[1], data2[2]))

    r1_out.communicate()
    r2_out.communicate()
    r1.close()
    r2.close()


def add_umi_se(input_dir, output_dir):
    R1 = glob.glob(input_dir + '/*R1*fastq*')[0]
    umi = glob.glob(input_dir + '/*R2*fastq*')[0]
    sample_name = R1.split('/')[-1].split('_')[0]
    r1 = open(output_dir + '/' + sample_name + '_R1_UMI_001.fastq.gz', 'w')
    # Use subprocess to pipe to gzip
    r1_out = subprocess.Popen('gzip', stdin=subprocess.PIPE, stdout=r1)

    with io.BufferedReader(gzip.open(R1, 'rb')) as read1, \
            io.BufferedReader(gzip.open(umi, 'rb')) as uid:
        for data1, umidata in itertools.izip(readfq(read1),
                                             readfq(uid)):
            r1_name, r1_suffix = data1[0].split()
            umi_name = umidata[0].split()[0]
            umi_id = umidata[1]
            assert(r1_name == umi_name), 'Mismatch in read names'
            r1_umi = r1_name + ':' + umi_id + ' ' + r1_suffix
            writefq(r1_out.stdin, (r1_umi, data1[1], data1[2]))

    r1_out.communicate()
    r1.close()

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', help='Directory with Fastqs')
    parser.add_argument('-o', '--output_dir', help='Output directory')
    parser.add_argument('-s', '--single_end', action='store_true',
                        default=False)

    opts = parser.parse_args()
    if opts.single_end:
        add_umi_se(opts.input_dir, opts.output_dir)
    else:
        add_umi(opts.input_dir, opts.output_dir)
