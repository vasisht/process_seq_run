#!/usr/bin/python

import sys
import gzip
from collections import defaultdict
import argparse


def count_barcodes(infile, barcodes):

    inf = gzip.open(infile, 'r')

    header = inf.readline()

    while header:
        seq = inf.readline()
        comment = inf.readline()
        quality = inf.readline()
        barcode = header.rstrip().split(':')[-1]
        barcodes[barcode] += 1

        header = inf.readline()

    inf.close()
    return barcodes


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs='+',
                        help="Undetermined fastq files")
    parser.add_argument("-c", "--count", type=int, default=15)

    opts = parser.parse_args()

    barcodes = defaultdict(int)
    for f in opts.infiles:
        barcodes = count_barcodes(f, barcodes)

    for k in sorted(barcodes.keys(), key=barcodes.get, reverse=True)[:opts.count]:
        print k, barcodes[k]
