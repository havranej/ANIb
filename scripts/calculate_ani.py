import pandas as pd
import os

INPUTS = snakemake.input
OUTPUT_PATH = snakemake.output[0]

results = []

for input in INPUTS:
    blast_results = pd.read_csv(input, sep = "\t")
    ani = blast_results.aln_percent_identity.mean()
    
    basename = os.path.basename(input)
    basename = basename.strip(".blast.tsv")
    genomes = basename.split("_vs_")

    results.append({"Seq1": genomes[0], "Seq2": genomes[1], "ANI": ani})

ani_table = pd.DataFrame(results)
ani_table.to_csv(OUTPUT_PATH)
