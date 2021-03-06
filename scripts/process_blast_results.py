
import pandas as pd

INPUT_PATH = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]

EVAL_THR = snakemake.params["eval_threshold"]
COVERAGE_THR = snakemake.params["coverage_threshold"]
IDENTITY_THR = snakemake.params["identity_threshold"]

blast_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None, names = ["qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "aln_length", "gaps", "aln_percent_identity", "nident", "evalue"])
blast_best = blast_results.copy().iloc[blast_results.groupby('qseqid').nident.idxmax()]

blast_best["aligned_pairs"] = blast_best.aln_length - blast_best.gaps
blast_best["query_coverage"] = blast_best.aligned_pairs / blast_best.qlen
blast_best["query_percent_identity"] = blast_best.nident / blast_best.qlen

# Filtering as in Goris et al. 2007
blast_best_filtered = blast_best[(blast_best.query_coverage > COVERAGE_THR) & (blast_best.query_percent_identity > IDENTITY_THR) & (blast_best.evalue < EVAL_THR)]

blast_best_filtered.to_csv(OUTPUT_PATH, sep = "\t")
