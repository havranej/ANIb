# ANIb

A Snakemake pipeline for ANI computation with the use of BLAST

## Parameters
* `input_dir`: A path to a directory with genomes. All vs. all computation will be performed.
* `output_dir`
* `chunk_size`: Length of query genome fragments which are then searched against the reference
* `blast_threads`: Number of threads for an individual BLAST search
* `blast_eval`: E-value filter in BLAST
