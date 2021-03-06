import gzip
from pathlib import Path


import pytest
import matcha

from .utils import hamming_dist, random_sequence

def run_matcher(f):
    """Setup and run a matcher for the fastq data
    
    Args:
        f (FastqReader): FastqReader handler with the relevant reads
    """
    # Make it so reads 1 and 4 match. Read 1 with 1 mismatch in the barcode 
    i5_matcher = matcha.ListMatcher(["TCCGAGCC", "ACAGGCGC"], ["i5_1", "i5_4"])
    i7_matcher = matcha.ListMatcher(["GCCAATTC", "CTGTATTA"], ["i7_1", "i7_4"])

    f.add_barcode("cell_i5", i5_matcher, "I2")
    f.add_barcode("cell_i7", i7_matcher, "I1")

    f.set_output_names("{cell_i5}+{cell_i7}:{read_name}")
    while f.read_chunk(2):
        matches = (f.matches["cell_i5"].dist <= 1) & (f.matches["cell_i7"].dist <= 1)
        f.write_chunk(matches)
    f.close()


def test_multi_file_io(tmpdir):
    tmpdir = Path(str(tmpdir))
    f = matcha.FastqReader()
    for read in ["R1", "R2", "I1", "I2"]:
        path = tmpdir / read
        path.write_text(test_data[read])
        out_path = tmpdir / (read + "_out")
        f.add_sequence(read, path, out_path)

    run_matcher(f)

    for read in ["R1", "R2", "I1", "I2"]:
        output = (tmpdir / (read + "_out")).read_text().splitlines()
        input_text = test_data[read].splitlines()

        # Check that sequence and quality are preserved
        assert output[1:4] == input_text[1:4]
        assert output[5:8] == input_text[13:16]

        assert output[0] == f"@i5_1+i7_1:{input_text[0][1:]}"
        assert output[4] == f"@i5_4+i7_4:{input_text[12][1:]}"

def test_threaded_io(tmpdir):
    tmpdir = Path(str(tmpdir))
    f = matcha.FastqReader(threads=4)
    for read in ["R1", "R2", "I1", "I2"]:
        path = tmpdir / read
        path.write_text(test_data[read])
        out_path = tmpdir / (read + "_out")
        f.add_sequence(read, path, out_path)

    run_matcher(f)

    for read in ["R1", "R2", "I1", "I2"]:
        output = (tmpdir / (read + "_out")).read_text().splitlines()
        input_text = test_data[read].splitlines()

        # Check that sequence and quality are preserved
        assert output[1:4] == input_text[1:4]
        assert output[5:8] == input_text[13:16]

        assert output[0] == f"@i5_1+i7_1:{input_text[0][1:]}"
        assert output[4] == f"@i5_4+i7_4:{input_text[12][1:]}"

def test_gzip_io(tmpdir):
    tmpdir = Path(str(tmpdir))
    f = matcha.FastqReader()
    for read in ["R1", "R2", "I1", "I2"]:
        path = tmpdir / read
        gzip.open(path, "wt").write(test_data[read])
        out_path = tmpdir / (read + "_out.gz")
        f.add_sequence(read, path, out_path)

    run_matcher(f)

    for read in ["R1", "R2", "I1", "I2"]:
        output = gzip.open(tmpdir / (read + "_out.gz"), 'rt').read().splitlines()
        input_text = test_data[read].splitlines()

        # Check that sequence and quality are preserved
        assert output[1:4] == input_text[1:4]
        assert output[5:8] == input_text[13:16]

        assert output[0] == f"@i5_1+i7_1:{input_text[0][1:]}"
        assert output[4] == f"@i5_4+i7_4:{input_text[12][1:]}"

test_data = {}
test_data["I1"] = """\
@NB551514:265:H5KHFBGXC:1:23208:10434:9061 1:N:0:0
GCCAATTC
+
AAAAAEEE
@NB551514:265:H5KHFBGXC:1:12106:23211:12984 1:N:0:0
CGTACTAG
+
AAAAAEEE
@NB551514:265:H5KHFBGXC:1:23207:19364:4497 1:N:0:0
CTCATGGG
+
A/AAA/AE
@NB551514:265:H5KHFBGXC:1:21112:8047:14790 1:N:0:0
CTGTATTA
+
AAAAAEEE
@NB551514:265:H5KHFBGXC:1:21105:9516:13053 1:N:0:0
ATCACTCG
+
AAA/AAEA
"""

test_data["I2"] = """\
@NB551514:265:H5KHFBGXC:1:23208:10434:9061 2:N:0:0
TCCGTGCC
+
AAAAAEEE
@NB551514:265:H5KHFBGXC:1:12106:23211:12984 2:N:0:0
GCGATCTA
+
AAAAAEEE
@NB551514:265:H5KHFBGXC:1:23207:19364:4497 2:N:0:0
ATCATGTT
+
A//AA/EA
@NB551514:265:H5KHFBGXC:1:21112:8047:14790 2:N:0:0
ACAGGCGC
+
6A6AA6EE
@NB551514:265:H5KHFBGXC:1:21105:9516:13053 2:N:0:0
TGTAGATT
+
AAAAA6EE
"""

test_data["R1"] = """\
@NB551514:265:H5KHFBGXC:1:23208:10434:9061 1:N:0:0
TCATTTGCGTGCCGAGTAAAATGTCCGCTTTTCTGT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB551514:265:H5KHFBGXC:1:12106:23211:12984 1:N:0:0
GGTCATGAAGGCCACCTATCCCAAGTGAAATTCTGA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB551514:265:H5KHFBGXC:1:23207:19364:4497 1:N:0:0
TCGTCACGGACGCCAATCAGCAGGAATACCGATCGA
+
A/AAAE6//E/A///EE6E//E//E/E///AEE66E
@NB551514:265:H5KHFBGXC:1:21112:8047:14790 1:N:0:0
AAGGATGATTTTTTTTTTTTTTTTTTTTTTTTTTTC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEE//
@NB551514:265:H5KHFBGXC:1:21105:9516:13053 1:N:0:0
CACGAGGACATCGACGCCGACACGATCAACGCGGTG
+
AAAAA/EEEEEEAEEEAEEEEE/EEEEEEEAEAEEE
"""

test_data["R2"] = """\
@NB551514:265:H5KHFBGXC:1:23208:10434:9061 2:N:0:0
TAAACGAGTTTGGCGACAGAAAAGCGGACATTTTAC
+
AAA6AEEEEEEEEEEEEEEEEAAEEEEEE6EEEEEE
@NB551514:265:H5KHFBGXC:1:12106:23211:12984 2:N:0:0
ATCTCATACCATCACCTTTGGATGAAGGGTCATCAG
+
AAAAAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEE
@NB551514:265:H5KHFBGXC:1:23207:19364:4497 2:N:0:0
ACGTCGAAAGGATGCTGGTTCGATCTGGAGTCATGC
+
AAA///EA/E/EE6/////6//EA////E/E/////
@NB551514:265:H5KHFBGXC:1:21112:8047:14790 2:N:0:0
AAAGTCACTCTGCCGGAAAAAAAAAAAAAAAAAAAA
+
AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEE
@NB551514:265:H5KHFBGXC:1:21105:9516:13053 2:N:0:0
TCCTCGAGCACCGCGTTGATCGTGTCGGCGTCGATG
+
AAAAAEEE//AEEEEEEEEEEEEEEEEEEEEEE6EE
"""
