#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""
This script converts two types of constructs to a unified format
SYNTAX: convert_to_cagi.py <trivariate_file> <fasta_reference_1>, <fasta_reference_2>, ...
FASTA IDs must be "GENE_ID::chrN:start-end"
Output is saved in the same directory, names: "{source_name}.GENE_ID.unified.tsv"
"""

import numpy as np
import pandas as pd
import os
import re
from sys import argv
from Bio import SeqIO

INPUT_FILE_NAME = argv[1]
FASTA = argv[2:]
# BED = argv[3]

FIELDS = ["chr", "pos", "ref", "alt", "val", "conf"]
FIELDS_TRIV = ["name", "pos", "alt", "val", "conf"]


sheet = pd.read_csv(INPUT_FILE_NAME, comment="#", sep="\t", header=None, names=FIELDS_TRIV)

reference = dict()
for file in FASTA:
    for record in SeqIO.parse(file, "fasta"):
        # "ECR11::chr2:169939081-169939701"
        name = record.id.split("::")[0].replace("_", "-")
        start = int(re.search(r"(?<=:)(\d+)(?=-)", record.id).group(0))
        # 0-based; becomes 1-based with summation
        chrom = re.search(r"(?<=chr)[XY\d]{1,2}(?=:)", record.id).group(0)
        reference[name] = (chrom, start, str(record.seq))

for id in sorted(set(sheet["name"])):
    out_sheet = pd.DataFrame(columns=FIELDS)
    id_sheet = sheet[sheet["name"] == id]
    chrom, start, refseq = reference[id.upper()]
    for row_n, old_row in id_sheet.iterrows():
        row = pd.Series(index=FIELDS, name=row_n, dtype=object)
        n = int(old_row["pos"])
        row["chr"] = chrom
        row["pos"] = start + n
        row["ref"] = refseq[n - 1]
        row["alt"] = old_row["alt"]
        row["val"] = old_row["val"]
        row["conf"] = old_row["conf"]
        out_sheet.loc[row_n] = row
        if row["ref"] == row["alt"]:
            print(row)
            raise ValueError("reference must not be equal to alternative")
    outfile_name = "{}.{}.unified.tsv".format(os.path.splitext(INPUT_FILE_NAME)[0], id.upper())
    out_sheet.to_csv(outfile_name, sep='\t', index=False)

print("done")