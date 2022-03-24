#!/usr/bin/python
import sys
import argparse
from collections import defaultdict
from rank_barcodes import count_barcodes
import distance

i5_barcodes = {'A501':['TGAACCTT','AAGGTTCA'],
               'A502':['TGCTAAGT','ACTTAGCA'],
               'A503':['TGTTCTCT','AGAGAACA'],
               'A504':['TAAGACAC','GTGTCTTA'],
               'A505':['CTAATCGA','TCGATTAG'],
               'A506':['CTAGAACA','TGTTCTAG'],
               'A507':['TAAGTTCC','GGAACTTA'],
               'A508':['TAGACCTA','TAGGTCTA'],
               'D501':['TATAGCCT', 'AGGCTATA'],
               'D502':['ATAGAGGC', 'GCCTCTAT'],
               'D503':['CCTATCCT', 'AGGATAGG'],
               'D504':['GGCTCTGA', 'TCAGAGCC'],
               'D505':['AGGCGAAG', 'CTTCGCCT'],
               'D506':['TAATCTTA', 'TAAGATTA'],
               'D507':['CAGGACGT', 'ACGTCCTG'],
               'D508':['GTACTGAC', 'GTCAGTAC']}

i7_barcodes = {'A701':'ATCACGAC',
               'A702':'ACAGTGGT',
               'A703':'CAGATCCA',
               'A704':'ACAAACGG',
               'A705':'ACCCAGCA',
               'A706':'AACCCCTC',
               'A707':'CCCAACCT',
               'A708':'CACCACAC',
               'A709':'GAAACCCA',
               'A710':'TGTGACCA',
               'A711':'AGGGTCAA',
               'A712':'AGGAGTGG',
               'D701':'ATTACTCG',
               'D702':'TCCGGAGA',
               'D703':'CGCTCATT',
               'D704':'GAGATTCC',
               'D705':'ATTCAGAA',
               'D706':'GAATTCGT',
               'D707':'CTGAAGCT',
               'D708':'TAATGCGC',
               'D709':'CGGCTATG',
               'D710':'TCCGCGAA',
               'D711':'TCTCGCGC',
               'D712':'AGCGATAG'}

i5_list = i5_barcodes.values()
i5 = [i for s in i5_list for i in s]
i7 = i7_barcodes.values()

def get_distance(barcode, barcode_list):
    dist = []
    for b in barcode_list:
        d = distance.hamming(barcode, b)
        # Ignore the distance from itself
        if d != 0:
            dist.append(d)
    return dist

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--summary",
                    help="index summary from Run")
parser.add_argument("-f", "--fastq", nargs='+',
                    help="Undetermined fastq files")

opts = parser.parse_args()

sys.stderr.write("Parsing undetermined barcodes\n")
undetermined_barcodes = defaultdict(int)
for f in opts.fastq:
    undetermined_barcodes = count_barcodes(f, undetermined_barcodes)

sys.stderr.write("Parsing index summary file\n")
mapped_barcodes = defaultdict(int)
with open(opts.summary, "r") as inf:
    for lines in inf:
        if lines[0] == "#" or lines[0] == "L":
            continue
        else:
            fields = lines.rstrip().split(",")
            if fields[2] == '1':
                mapped_barcodes[fields[3]] += int(fields[6])

all_barcodes = {}
all_barcodes.update(mapped_barcodes)
all_barcodes.update(undetermined_barcodes)
total_reads = sum(all_barcodes.values())

barcode_perc = {}
for k,v in all_barcodes.iteritems():
    barcode_perc[k] = float(v*100.0)/float(total_reads)

filtered_barcodes = {k: v for k,v in barcode_perc.items() if v >= 0.05}

sys.stderr.write("Number of filtered barcodes %d\n" % len(filtered_barcodes.keys()))

barcode_dist = {}
# Get hamming distance to barcodes in samplesheet
for b in filtered_barcodes.keys():
    barcode_dist[b] = get_distance(b, mapped_barcodes.keys())

print "Barcode\t%Reads\tMinDist\tSampleSheet\tIlluminaI7\tIlluminaI5"

for k in filtered_barcodes.keys():
    is_i5, is_i7 = False, False
    k_i7, k_i5 = k.split("+")
    if k_i7 in i7:
        is_i7 = True
    if k_i5 in i5:
        is_i5 = True
    print "%s\t%.2f\t%d\t%s\t%s\t%s" % (k, barcode_perc[k], min(barcode_dist[k]),
                                mapped_barcodes.has_key(k), is_i7, is_i5)
