#!/usr/bin/env python
import gzip

import pandas as pd
import numpy as np
import matcha
import tqdm

# Read layout -- 
# R1: 16bp cell barcode, then 12bp UMI
# R2: 10bp of spacer, then ADT barcode
# I1: 8bp sample barcode

# Output for each read 
# - Sample barcode name + mismatch count
# - UMI sequence
# - Cell barcode best match, dist, and second_best_dist
# - ADT barcode best match, dist, and second_best_dist

output = open("pbmc_1k_protein_v3/example_matching_stats.csv", "w", newline="")

f = matcha.FastqReader(threads=4)
f.add_sequence("R1", "pbmc_1k_protein_v3/R1.fastq.gz")
f.add_sequence("R2", "pbmc_1k_protein_v3/R2.fastq.gz")
f.add_sequence("I1", "pbmc_1k_protein_v3/I1.fastq.gz")

sample_barcode = matcha.ListMatcher(
    ["CAGTACTG", "AGTAGTCT", "GCAGTAGA", "TTCCCGAC"],
    ["SI-GA-A3", "SI-GA-A3", "SI-GA-A3", "SI-GA-A3"]
)
f.add_barcode("sample", sample_barcode, "I1")

valid_10x_cell_barcodes = pd.read_csv(
    "pbmc_1k_protein_v3/3M-february-2018.txt.gz",
    sep="\t",
    names=["gene_barcode", "feature_barcode"]
)
cell_barcode = matcha.HashMatcher(
    valid_10x_cell_barcodes["feature_barcode"],
    max_mismatches=1,
    subsequence_count=2
)
f.add_barcode("cell", cell_barcode, "R1")

adt_barcode_list = pd.read_csv(
    "pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv"
)
adt_barcode = matcha.ListMatcher(
    adt_barcode_list["sequence"],
    adt_barcode_list["id"]
)
f.add_barcode("adt", adt_barcode, "R2", match_start=10)

first_write = True
total = None
chunk_size = 10000
progress = tqdm.tqdm(disable=None)
while f.read_chunk(chunk_size):
    pass_filter = (f.get_match_result("sample", "dist") <= 1) & \
        (f.get_match_result("sample", "second_best_dist") > 1)

    df = pd.DataFrame({
        "cell": f.get_match_result("cell", "label"),
        "UMI": f.get_sequence_read("R1", start=16),
        "ADT": f.get_match_result("adt", "label"),
        "cell_dist": f.get_match_result("cell", "dist"),
        "cell_second_best_dist": f.get_match_result("cell", "second_best_dist"),
        "adt_dist": f.get_match_result("adt", "dist"),
        "adt_second_best_dist": f.get_match_result("adt", "second_best_dist"),
    })

    df[pass_filter].to_csv(
        output,
        header=first_write,
        index=False)

    progress.update(chunk_size)
    first_write=False