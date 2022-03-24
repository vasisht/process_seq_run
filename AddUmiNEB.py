#!/usr/bin/python
import gzip
import glob
import itertools
import io
import subprocess
from readfq import readfq

'''
Commands to generate read and barcode fastq files from NEBNext Direct
java -Xmx4G -jar /usr/local/bin/picard.jar ExtractIlluminaBarcodes 
BASECALLS_DIR=Data/Intensities/BaseCalls/ LANE=001 READ_STRUCTURE=300T8B12M300T 
BARCODE_FILE=run_barcodes.txt METRICS_FILE=metrics_output.txt NUM_PROCESSORS=0

java -Xmx4G -jar /usr/local/bin/picard.jar IlluminaBasecallsToFastq  READ_STRUCTURE=300T8B12M300T 
BASECALLS_DIR=Data/Intensities/BaseCalls/ LANE=001 MULTIPLEX_PARAMS=barcodes.txt NUM_PROCESSORS=6 
COMPRESS_OUTPUTS=TRUE READ_NAME_FORMAT=ILLUMINA RUN_BARCODE=BM92W MACHINE_NAME=M03295 
INCLUDE_NON_PF_READS=FALSE IGNORE_UNEXPECTED_BARCODES=TRUE
'''


def writefq(fh, data):
    fh.write('@' + data[0] + '\n')
    fh.write(data[1] + '\n')
    fh.write('+\n')
    fh.write(data[2] + '\n')


def add_umi_neb(input_dir, output_dir):
    R1 = glob.glob(input_dir + '/*.1.fastq.gz')[0]
    R2 = glob.glob(input_dir + '/*.2.fastq.gz')[0]
    umi = glob.glob(input_dir + '/*index*')[0]
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

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', help='Directory with Fastqs')
    parser.add_argument('-o', '--output_dir', help='Output directory')

    opts = parser.parse_args()
    add_umi_neb(opts.input_dir, opts.output_dir)
