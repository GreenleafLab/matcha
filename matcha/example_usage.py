#Example use of Matcher library

cell_i5 = Matcher(sequences, labels)
cell_i7 = Matcher(sequences, labels)

f = FastqReader()
f.add_sequence("R1", input_path="in/R1.fastq.gz", output_path="out/R1.fastq.gz")
f.add_sequence("R2", input_path="in/R2.fastq.gz", output_path="out/R2.fastq.gz")

f.add_barcode("cell_i5", cell_i5, "R1")
f.add_barcode("cell_i7", cell_i7, "R2")

f.set_output_names("{cell_i5}+{cell_i7}:{lane}:{tile}:{x}:{y}")

while f.read_chunk():
    pass_filter = (f.matches.cell_i5.dist <= 1) & (f.matches.cell_i7.dist <= 1)
    f.write_chunk(filter=pass_filter)
