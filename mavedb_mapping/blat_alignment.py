from Bio import SearchIO
import pandas as pd
import subprocess

from mavedb_mapping import path_to_hg38_file

"""
    Module for performing BLAT Alignment using target sequence.
    Returns Query sequence locations with corresponding locations in the
    hit sequence obtained by BLAT.

"""


def extract_blat_output(dat: dict):
    """
    Runs a BLAT Query and returns the output.

    Parameters
    ----------
        dat: dict
            Dictionary containing data required for mapping.

    Returns:
    --------
        BLAT Output
    """
    blat_file = open("blat_query.fa", "w")
    blat_file.write(">" + "query" + "\n")
    blat_file.write(dat["target_sequence"] + "\n")
    blat_file.close()
    # minimum match 50%
    min_score = len(dat["target_sequence"]) // 2

    if dat["target_sequence_type"] == "protein":
        command = f"blat {path_to_hg38_file} -q=prot -t=dnax -minScore={min_score} blat_query.fa blat_out.psl"
        process = subprocess.run(command, shell=True)
    else:
        command = (
            f"blat {path_to_hg38_file} -minScore={min_score} blat_query.fa blat_out.psl"
        )
        process = subprocess.run(command, shell=True)
    try:
        output = SearchIO.read("blat_out.psl", "blat-psl")
    except ValueError:
        try:
            command = f"blat {path_to_hg38_file} -q=dnax -t=dnax -minScore={min_score} blat_query.fa blat_out.psl"
            process = subprocess.run(command, shell=True)
            output = SearchIO.read("blat_out.psl", "blat-psl")
        except ValueError:
            return None
    return output


def obtain_hit_starts(output, hit: int):
    # a hit is a portion of similarity between query seq and matched seq
    """
    Helper function to obtain HSP.
    Returns the starts of hit sequence.
    """
    hit_starts = list()
    for n in range(len(output[hit])):
        hit_starts.append(output[hit][n].hit_start)
    return hit_starts


def obtain_hsp(output):
    """
    Obtains high-scoring pairs (HSP) for query sequence.

    Parameters
    ----------
        output: BLAT output

    Returns
    -------
        HSP

    """
    hit_dict = {}
    for c in range(len(output)):
        max_hit = output[c][0].score
        for e in range(len(output[c])):
            if output[c][e].score > max_hit:
                max_hit = output[c][e].score
        hit_dict[c] = max_hit

    # Using top scoring hit
    hit = max(hit_dict, key=hit_dict.get)
    hit_starts = obtain_hit_starts(output, hit)
    hsp = output[hit][hit_starts.index(max(hit_starts))]
    return hsp


def get_query_and_hit_ranges(hsp, output) -> tuple:
    """

    Parameters
    ----------
        hsp:
            High scoring pairs obtained by using top scoring hit
        output:
            BLAT output
        gsymb: str
            Gene symbol

    Returns
    -------
        Tuple containing the chromosome, strand, coverage, query ranges, and hit ranges.

    """
    strand = hsp[0].query_strand
    coverage = 100 * (hsp.query_end - hsp.query_start) / output.seq_len
    query_ranges = []
    hit_ranges = []
    for j in range(len(hsp)):
        test_file = open("blat_output_test.txt", "w")
        test_file.write(str(hsp[j]))
        test_file.close()

        query_string = ""
        hit_string = ""

        test_file = open("blat_output_test.txt", "r")
        for k, line in enumerate(test_file):
            if k == 1:
                chrom = line.strip("\n")
            if k == 2:
                query_string = line.strip("\n")
            if k == 3:
                hit_string = line.strip("\n")
        test_file.close()

        chrom = chrom.split(" ")[9].strip("chr")
        query_string = query_string.split(" ")
        hit_string = hit_string.split(" ")
        query_ranges.append(query_string[2])
        hit_ranges.append(hit_string[4])

    return chrom, strand, coverage, query_ranges, hit_ranges


def mave_to_blat(dat: dict) -> dict:
    """

    Performs BLAT Alignment on MaveDB scoreset data.

    Parameters
    ----------
        dat: dict
            Dictionary containing data from MaveDB scoresets.

    Returns
    -------
        mave_blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment

    """
    output = extract_blat_output(dat)
    if output is not None:
        hsp = obtain_hsp(output)
        (
            chrom,
            strand,
            coverage,
            query_ranges,
            hit_ranges,
        ) = get_query_and_hit_ranges(hsp, output)
        qh_dat = {"query_ranges": query_ranges, "hit_ranges": hit_ranges}
        qh_dat = pd.DataFrame(data=qh_dat)
        mave_blat_dict = {
            "chrom": chrom,
            "strand": strand,
            "target_type": dat["target_type"],
            "coverage": coverage,
            "hits": qh_dat,
        }

    else:
        qh_dat = {"query_ranges": ["NA"], "hit_ranges": ["NA"]}
        qh_dat = pd.DataFrame(data=qh_dat)
        mave_blat_dict = {
            "chrom": "NA",
            "strand": "NA",
            "target_type": "NA",
            "coverage": "NA",
            "hits": qh_dat,
        }

    return mave_blat_dict
