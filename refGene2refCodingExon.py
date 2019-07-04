#!/usr/bin/env python3

import argparse
import pandas as pd
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("-r", dest="refgene", required=True)
args = parser.parse_args()

# canonical chromosomes
canonical = [str(i) for i in range(1, 23)] + ["X", "Y", "M"]


def main(refgene):
    with open("refCodingExon.unsorted.bed", "w") as fo:
        header = (
            "chr\tstart\tend\tinfo\tstrand\tprev_exon_start|end\tnext_exon_start|end\n"
        )
        fo.write(header)
        d = isoform_cdslen_dict(refgene)
        with open(refgene) as fi:
            for line in fi:
                try:
                    fo.write(parse_refgene_line(line, d) + "\n")
                except:
                    pass

    df = pd.read_csv("refCodingExon.unsorted.bed", sep="\t")

    # sort by chr and start
    df.sort_values("chr", inplace=True)
    df = df.groupby("chr").apply(pd.DataFrame.sort_values, "start")

    df["chr"] = df.apply(format_chromosom, axis=1)

    df.to_csv("refCodingExon.bed", sep="\t", index=False, header=False)

    sp.call(["rm", "refCodingExon.unsorted.bed"])


def isoform_cdslen_dict(refgene):
    """Makes a dict {RefSeq accession: len of cds}
    """
    d = {}
    with open(refgene) as f:
        for line in f:
            lst = line.rstrip().split("\t")
            acc = lst[1]
            chr = lst[2]
            strand = lst[3]

            # convert to 1-based coordinate for start
            cds_start = int(lst[6]) + 1
            exon_starts = [int(pos) + 1 for pos in lst[9].rstrip(",").split(",")]

            # no need to convert for end
            cds_end = int(lst[7])
            exon_ends = [int(pos) for pos in lst[10].rstrip(",").split(",")]

            # frame info to identify coding exons (-1 for UTR exons)
            frames = [int(frame) for frame in lst[15].rstrip(",").split(",")]

            # count the number of coding exons
            total_num_of_exons = len(lst[15].replace("-1,", "").rstrip(",").split(","))

            # reverse positions for negative stand starts/ends
            if strand == "-":
                exon_starts = exon_starts[::-1]
                exon_ends = exon_ends[::-1]
                frames = frames[::-1]

            cds_len = 0
            this_is_i_th_coding_exon = 1
            for start, end, frame in zip(exon_starts, exon_ends, frames):
                # skip if the entire exon is UTR
                if frame == -1:
                    pass
                # coding exons on positive strand
                elif frame != -1 and strand == "+":
                    # start
                    if this_is_i_th_coding_exon == 1:
                        coding_exon_i_start = cds_start
                    else:
                        coding_exon_i_start = start
                    # end
                    if this_is_i_th_coding_exon == total_num_of_exons:
                        coding_exon_i_end = cds_end
                    else:
                        coding_exon_i_end = end

                    # cds spans
                    if this_is_i_th_coding_exon == 1:
                        span = int(end) - int(cds_start) + 1
                    elif this_is_i_th_coding_exon == total_num_of_exons:
                        span = int(cds_end) - int(start) + 1
                    else:
                        span = int(end) - int(start) + 1

                    cds_len += span
                    this_is_i_th_coding_exon += 1

                # coding exons on negative strand
                elif frame != -1 and strand == "-":
                    # start
                    if this_is_i_th_coding_exon == total_num_of_exons:
                        coding_exon_i_start = cds_start
                    else:
                        coding_exon_i_start = start
                    # end
                    if this_is_i_th_coding_exon == 1:
                        coding_exon_i_end = cds_end
                    else:
                        coding_exon_i_end = end

                    # cds spans
                    if this_is_i_th_coding_exon == total_num_of_exons:
                        span = int(end) - int(cds_start) + 1
                    elif this_is_i_th_coding_exon == 1:
                        span = int(cds_end) - int(start) + 1
                    else:
                        span = int(end) - int(start) + 1

                    cds_len += span
                    this_is_i_th_coding_exon += 1

            # record results for isoforms on canonical chromosomes
            if chr.lstrip("chr") in canonical:
                d[acc] = cds_len
            else:
                pass

    return d


def format_chromosom(row):
    if row["chr"] == 23:
        return "chrX"
    elif row["chr"] == 24:
        return "chrY"
    elif row["chr"] == 25:
        return "chrM"
    else:
        return "chr" + str(row["chr"])


def parse_refgene_line(line, iso_len_dict):
    lst = line.rstrip().split("\t")
    acc = lst[1]
    chromosome = lst[2]

    # remove non-canonicals and format chromosomes
    chr = chromosome.replace("chr", "")
    if chr not in canonical:
        return None
    if chr == "X":
        chr = "23"
    elif chr == "Y":
        chr = "24"
    elif chr == "M":
        chr = "25"

    strand = lst[3]

    # convert to 1-based coordinate for start
    cds_start = int(lst[6]) + 1
    exon_starts = [int(pos) + 1 for pos in lst[9].rstrip(",").split(",")]

    # no need to convert for end
    cds_end = int(lst[7])
    exon_ends = [int(pos) for pos in lst[10].rstrip(",").split(",")]

    # frame info to identify coding exons (-1 for UTR exons)
    frames = [int(frame) for frame in lst[15].rstrip(",").split(",")]

    gene_symbol = lst[12]

    # count the number of coding exons
    total_num_of_exons = len(lst[15].replace("-1,", "").rstrip(",").split(","))

    # reverse positions for negative stand starts/ends
    if strand == "-":
        exon_starts = exon_starts[::-1]
        exon_ends = exon_ends[::-1]
        frames = frames[::-1]

    parsed = []
    from_start_codon = 1
    this_is_i_th_coding_exon = 1
    for start, end, frame in zip(exon_starts, exon_ends, frames):
        # skip if the entire exon is UTR
        if frame == -1:
            pass
        # coding exons
        else:
            # on positive strand
            if strand == "+":
                # start
                if this_is_i_th_coding_exon == 1:
                    coding_exon_i_start = cds_start
                else:
                    coding_exon_i_start = start
                # end
                if this_is_i_th_coding_exon == total_num_of_exons:
                    coding_exon_i_end = cds_end
                else:
                    coding_exon_i_end = end
                # cds spans
                if this_is_i_th_coding_exon == 1:
                    span = int(end) - int(cds_start) + 1
                elif this_is_i_th_coding_exon == total_num_of_exons:
                    span = int(cds_end) - int(start) + 1
                else:
                    span = int(end) - int(start) + 1

            # on negative strand
            else:
                # start
                if this_is_i_th_coding_exon == total_num_of_exons:
                    coding_exon_i_start = cds_start
                else:
                    coding_exon_i_start = start
                # end
                if this_is_i_th_coding_exon == 1:
                    coding_exon_i_end = cds_end
                else:
                    coding_exon_i_end = end
                # cds spans
                if this_is_i_th_coding_exon == total_num_of_exons:
                    span = int(end) - int(cds_start) + 1
                elif this_is_i_th_coding_exon == 1:
                    span = int(cds_end) - int(start) + 1
                else:
                    span = int(end) - int(start) + 1

            cds_len_of_this_isoform = iso_len_dict[acc]
            ith_coding_exon = (
                chr
                + "\t"
                + str(coding_exon_i_start)
                + "\t"
                + str(coding_exon_i_end)
                + "\t"
                + acc
                + "|"
                + gene_symbol
                + "|"
                + str(this_is_i_th_coding_exon)
                + "|"
                + str(total_num_of_exons)
                + "|"
                + str(from_start_codon)
                + "|"
                + str(cds_len_of_this_isoform)
                + "\t"
                + strand
            )

            parsed.append(ith_coding_exon)

            this_is_i_th_coding_exon += 1
            from_start_codon += span

    # adding i-1 and i+1 exon starts
    if len(parsed) == 0:
        pass
    elif len(parsed) == 1:
        exon = parsed[0]
        data = exon + "\t-1|-1\t-1|-1"
        return data + "\n"
    else:
        parsed_added = []
        for i in range(len(parsed)):
            if i == 0:
                second_coding_exon = parsed[1]
                second_coding_exon_start = second_coding_exon.split("\t")[1]
                second_coding_exon_end = second_coding_exon.split("\t")[2]
                first_coding_exon = (
                    parsed[0]
                    + "\t"
                    + "-1|-1"
                    + "\t"
                    + second_coding_exon_start
                    + "|"
                    + second_coding_exon_end
                )
                parsed_added.append(first_coding_exon)
            elif 0 < i < len(parsed) - 1:
                prev_coding_exon = parsed[(i - 1)]
                prev_coding_exon_start = prev_coding_exon.split("\t")[1]
                prev_coding_exon_end = prev_coding_exon.split("\t")[2]
                next_coding_exon = parsed[(i + 1)]
                next_coding_exon_start = next_coding_exon.split("\t")[1]
                next_coding_exon_end = next_coding_exon.split("\t")[2]
                current_coding_exon = (
                    parsed[i]
                    + "\t"
                    + prev_coding_exon_start
                    + "|"
                    + prev_coding_exon_end
                    + "\t"
                    + next_coding_exon_start
                    + "|"
                    + next_coding_exon_end
                )
                parsed_added.append(current_coding_exon)
            else:
                second_last_coding_exon = parsed[-2]
                second_last_coding_exon_start = second_last_coding_exon.split("\t")[1]
                second_last_coding_exon_end = second_last_coding_exon.split("\t")[2]
                last_coding_exon = (
                    parsed[-1]
                    + "\t"
                    + second_last_coding_exon_start
                    + "|"
                    + second_last_coding_exon_end
                    + "\t"
                    + "-1|-1"
                )
                parsed_added.append(last_coding_exon)
        return "\n".join(parsed_added)


if __name__ == "__main__":
    main(args.refgene)
