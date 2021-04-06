import pandas as pd
import os


QUERY_LENGTHS_PATH = snakemake.input[0]
INPUTS = snakemake.input[1:]
OUTPUT_PATH = snakemake.output[0]

results = []
query_lengths = pd.read_csv(QUERY_LENGTHS_PATH).drop_duplicates().set_index("genome")

for input_path in INPUTS:
    blast_results = pd.read_csv(input_path, sep = "\t", index_col = 0)
    ani = blast_results.aln_percent_identity.mean()
    
    basename = os.path.basename(input_path)
    basename = basename[:-10] # Remove last 10 characters (.blast.tsv)
    genomes = basename.split("_vs_")

    query_length = query_lengths.at[genomes[0], "length"]
    query_coverage = blast_results.aligned_pairs.sum() / query_length * 100

    results.append({"Seq1": genomes[0], "Seq2": genomes[1], "ANI": ani, "Coverage" : query_coverage})

ani_table = pd.DataFrame(results)
ani_table.to_csv(OUTPUT_PATH)
