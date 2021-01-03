Quick Example
=========================

This will give a quick example of matching ADT data from 10x CITE-seq dataset.
We'll be using the public pbmc_1k_protein_v3 dataset from 10x genomics available `here <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_protein_v3>`_,
downsampled to 50k reads. The full script and input dataset are available in the example folder of the github.

For each read in the input files, we'll match the sample, cell, and ADT barcodes,
and output match information as well as the UMI sequence in a csv file.

Setup
------
The read layout for this dataset is as follows: 

* R1: 16bp cell barcode, then 12bp UMI
* R2: 10bp of spacer, then ADT barcode
* I1: 8bp sample barcode

We'll start with the required imports, these aren't too important

.. code-block:: python

    import gzip

    import pandas as pd
    import numpy as np
    import matcha
    import tqdm

Opening files
--------------
Next, we open our output csv file

.. code-block:: python

    output = open("pbmc_1k_protein_v3/example_matching_stats.csv", "w", newline="")

Then, we specify our input files. First we create a :class:`~matcha.FastqReader` object.
This is the object that will perform all the actual reading and barcode matching.
Although we won't use it in this tutorial, it can also write back out filtered fastq files,
which you can use for demultiplexing or tagging read names with corrected single cell barcodes.

We'll make our reader use 3 threads, which will be used for reading files and
barcode matching in parallel

.. code-block:: python

    f = matcha.FastqReader(threads=3)

Next we specify our fastq file paths and add them to our FastqReader

.. code-block:: python

    f.add_sequence("R1", "pbmc_1k_protein_v3/R1.fastq.gz")
    f.add_sequence("R2", "pbmc_1k_protein_v3/R2.fastq.gz")
    f.add_sequence("I1", "pbmc_1k_protein_v3/I1.fastq.gz")

Barcode Sequences
-------------------

Now that we've added our sequences to our FastqReader object, it's time to 
specify our barcode sequences and where to find them in each read. Matcha has
two types of matching classess: 

#.  :class:`~matcha.ListMatcher`, which stores each valid barcode in a list, then
    compares each read to each barcode and picks the closest. This matcher can detect
    or correct an unlimited number of errors, but it becomes very slow if it needs to
    check too many valid sequences. As a rule of thumb, this is the best choice if you
    have under 20 valid sequences, and you should probably not use it to match over 100 sequences
#.  :class:`~matcha.HashMatcher`, which stores valid barcodes in a series of hash tables to
    allow efficient lookups of partial barcodes. This matcher can efficiently match against
    a huge number of valid barcodes, but it becomes slow if you need to be able to 
    correct more than about 10% of the bases in each barcode.

Note that the Matcher objects only store the valid barcodes and the matching strategy.
We specify what part of the sequence they should match against when we add them to our FastqReader object.

The first barcode we'll match is the sample barcode. In this case it's overkill
because our fastq files are already demultiplexed, but we'll do it for completeness.
Since there are only 4 valid barcode values, we'll use a ListMatcher:

.. code-block:: python

    sample_barcode = matcha.ListMatcher(
        ["CAGTACTG", "AGTAGTCT", "GCAGTAGA", "TTCCCGAC"],
        ["SI-GA-A3", "SI-GA-A3", "SI-GA-A3", "SI-GA-A3"]
    )
    f.add_barcode("sample", sample_barcode, "I1")

Next, we'll add the 10x barcodes. There are 3 million valid cell barcodes, so
we'll need to use the HashMatcher. We'll use the recommended parameters for single
basepair correction, although often it is a good idea to use only exact matches in 
this data. First we read in the valid barcodes from a file, then construct the
HashMatcher object and add it to match the start of R1.

.. code-block:: python

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

Finally, we'll add the ADT barcodes. Since there are only a handful we'll use
a ListMatcher. We'll read the valid barcodes from a file, then add the Matcher
to compare against R2 after skipping the first 10bp (since this 10x experiment
used Feature Barcoding with Total-Seq B format ADTs)

.. code-block:: python

    adt_barcode_list = pd.read_csv(
        "pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv"
    )
    adt_barcode = matcha.ListMatcher(
        adt_barcode_list["sequence"],
        adt_barcode_list["id"]
    )
    f.add_barcode("adt", adt_barcode, "R2", match_start=10)


Running the Matching
---------------------

So far, everything we've done has just been setup to configure how to read input
files, and how to match barcodes. Next, we'll do the actual reading. Matcha 
operates on chunks of reads in order to use fast numpy/pandas operations to greatly
speed up our processing. We'll use a chunk_size of 10,000 reads, which is sufficient
to not bottleneck on our python for loop.

All we'll do in the loop is test which reads pass filter (in this case, just
requiring a good match on the sample barcode). Then we'll make a pandas DataFrame
with the match information, and all the reads passing filter to our csv file.

.. code-block:: python

    first_write = True
    total = None
    chunk_size = 10000
    progress = tqdm.tqdm(disable=None)
    while f.read_chunk(chunk_size):
        pass_filter = f.get_match_result("sample", "dist") <= 1

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

Full script
------------
The full script for this example is available in the example folder of the `github <https://github.com/GreenleafLab/matcha/tree/master/example>`_.
On the example dataset, takes about a minute for me to run, most of which is spent indexing
the 3M cell barcodes.