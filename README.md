# ANIb

A Snakemake pipeline for ANI and query coverage computation with the use of BLAST

## Parameters
* `input_dir`: A path to a directory with genomes. 
* `input_list`: A path to a list of pairs. If this option is specified, only specified pairs will be compared. Otherwise, all vs. all computation will be performed.
* `output_dir`
* `input_extension`: Extension of the files with genomes (default: `fna`)
* `chunk_size`: Length of query genome fragments which are then searched against the reference (default: 1020)
* `blast_threads`: Number of threads for an individual BLAST search 
* `eval_filter`: Maximal E-value to pass the filters (default: 1.0E-15)
* `coverage_filter`: Minimal query coverage in reference to pass the filters (default: 0.7)
* `identity_filter`: Minimal fraction of query nucleotides aligned to reference as identical (default: 0.3)
* `blast_params`: Additional parameters to be passed to BLAST. `key: value` in the config file will be formatted as `-key value` in the command. Default values are `-dust no -xdrop_gap 150`