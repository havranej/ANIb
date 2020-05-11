# ANIb

A Snakemake pipeline for ANI computation with the use of BLAST

## Parameters
* `input_dir`: A path to a directory with genomes. All vs. all computation will be performed.
* `output_dir`
* `chunk_size`: Length of query genome fragments which are then searched against the reference (default: 1020)
* `blast_threads`: Number of threads for an individual BLAST search 
* `eval_threshold`: Maximal E-value to pass the filters (default: 1.0E-15)
* `coverage_threshold`: Minimal query coverage in reference to pass the filters (default: 0.7)
* `identity_threshold`: Minimal fraction of query nucleotides aligned to reference as identical (default: 0.3)
