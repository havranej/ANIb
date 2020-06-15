
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

INPUT_PATH = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]
CHUNK_SIZE = snakemake.params.chunk_size 

def split_into_chunks(seq, chunk_length = 1020):
    '''A generator to divide a sequence into chunks of n units.'''
    while seq:
        yield seq[:chunk_length]
        seq = seq[chunk_length:]

output_records = []

for record in SeqIO.parse(INPUT_PATH, "fasta"):
    for i, chunk in enumerate(split_into_chunks(record.seq, CHUNK_SIZE)):
        chunk_record = SeqRecord(chunk, id = record.id + "." + str(i), description = "")
        output_records.append(chunk_record)

with open(OUTPUT_PATH, "w") as output_handle:
    SeqIO.write(output_records, output_handle, "fasta")
