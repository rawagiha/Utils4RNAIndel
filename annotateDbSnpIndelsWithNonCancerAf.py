#!/usr/bin/env python3

import sys
import gzip
import pysam
import argparse
from variant import Variant


def main():
    args = get_args()

    dbsnp_vcf, gnomad_vcf, output_vcf = args.dbsnp, args.gnomad, args.output_vcf
    gnomad, genome = pysam.TabixFile(gnomad_vcf), pysam.FastaFile(args.fasta)

    fo = open(output_vcf, "w")

    fo.write(make_new_header(dbsnp_vcf, gnomad_vcf))

    for dbsnp_line in open_vcf(dbsnp_vcf):
        if dbsnp_line.startswith("#"):
            pass
        else:
            fo.write(annotate_dbsnp_line_with_non_cancer_af(dbsnp_line, gnomad, genome))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dbsnp", required=True)
    parser.add_argument("-g", "--gnomad", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-o", "--output-vcf", required=True)
    args = parser.parse_args()
    return args


def open_vcf(vcf):
    if vcf.endswith(".vcf"):
        f = open(vcf)
    elif vcf.endswith(".vcf.gz"):
        f = gzip.open(vcf, "rt")
    else:
        print("expected file extention is .vcf or .vcf.gz")
        sys.exit(1)
    return f


def make_new_header(dbsnp_vcf, gnomad_vcf):
    f = open_vcf(dbsnp_vcf)

    header = [line.rstrip() for line in f if line.startswith("#")]
    header1, header2 = header[: -1], header[-1]
    non_cancer_af_line = "##INFO=<ID=non_cancer_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in the gnomAD_non_cancer subset\">"
    input_source = "##gnomAD_file_used_for_annotation=" + gnomad_vcf
    header1.extend([non_cancer_af_line, input_source, header2])

    f.close()

    return "\n".join(header1) + "\n"


def make_indel(chrom, pos, ref, alt, genome):
    try:
        return Variant(chrom, pos, ref, alt, genome)
    except:
        return None


def parse_vcf_line(line, genome, with_locus_info=True, for_gnomad=False):
    lst = line.rstrip().split("\t")
    chrom, pos = lst[0], int(lst[1])
    ref, alts = lst[3], lst[4].split(",")
    info = lst[7].split(";")

    if for_gnomad:
        non_cancer_af = [
            float(i.replace("non_cancer_AF=", ""))
            for i in info
            if "non_cancer_AF=" in i
        ]
        if non_cancer_af:
            val = non_cancer_af[0]
        else:
            val = -1
        indels = [
            {"indel": make_indel(chrom, pos, ref, alt, genome), "non_cancer_af": val}
            for alt in alts
            if make_indel(chrom, pos, ref, alt, genome)
        ]
    else:
        indels = [
            make_indel(chrom, pos, ref, alt, genome)
            for alt in alts
            if make_indel(chrom, pos, ref, alt, genome)
        ]

    if with_locus_info:
        return chrom, pos, indels
    else:
        return indels


def annotate_dbsnp_line_with_non_cancer_af(dbsnp_line, gnomad, genome):
    chrom, pos, dbsnp_indels = parse_vcf_line(dbsnp_line, genome)

    window = 50

    non_cancer_af = []
    for gnomad_record in gnomad.fetch(chrom, pos - window, pos + window):
        gnomad_indels = parse_vcf_line(
            str(gnomad_record), genome, with_locus_info=False, for_gnomad=True
        )
        for dbsnp_indel in dbsnp_indels:
            for gnomad_indel in gnomad_indels:
                if dbsnp_indel == gnomad_indel["indel"]:
                    non_cancer_af.append(gnomad_indel["non_cancer_af"])

    af_data = ",".join([str(af) for af in non_cancer_af]) if non_cancer_af else ""
    
    if af_data:
        return dbsnp_line.rstrip() + ";" + "non_cancer_AF="+af_data + "\n"
    else:
        return dbsnp_line

if __name__ == "__main__":
    main()
