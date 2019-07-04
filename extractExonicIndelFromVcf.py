#!/usr/bin/env python3

import os
import sys
import gzip
import pysam
import argparse


def main():
    args = get_args()

    vcf, data_dir, output_vcf = args.vcf, args.data_dir, args.output_vcf

    if not vcf.endswith(".vcf") and not vcf.endswith(".vcf.gz") :
        print("expected file extention is .vcf or .vcf.gz")
        sys.exit(1)

    refgene = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
    exons = pysam.TabixFile(refgene)

    if vcf.endswith(".vcf"):
        f = open(vcf)
    else:
        f = gzip.open(vcf, "rt")

    fo = open(output_vcf, "w")

    for line in f:
        if line.startswith("#"):
            fo.write(line)
        else:
            lst = line.split("\t")
            if contains_indel(lst) and is_exonic(lst, exons):
                fo.write(line)
            else:
                pass
    f.close()
    fo.close()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcf", required=True)
    parser.add_argument("-d", "--data-dir", required=True)
    parser.add_argument("-o", "--output-vcf", required=True)
    args = parser.parse_args()
    return args


def contains_indel(lst):
    ref, alts = lst[3], lst[4]
    return sum([len(alt) != len(ref) for alt in alts.split(",")]) > 0


def is_exonic(lst, exons):
    """ may or may not be coding
    """

    chrom = lst[0] if lst[0].startswith("chr") else "chr" + lst[0]
    pos = int(lst[1])

    try:
        candidate_genes = exons.fetch(chrom, pos - 50, pos + 50)
    except:
        return False
    genes = [i for i in candidate_genes]
    
    return (genes != [])


if __name__ == "__main__":
    main()
