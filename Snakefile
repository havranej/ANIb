
import itertools
import os

from helpers import format_blast_params

OUT_DIR = config["output_dir"]
CHUNK_SIZE = config.get("chunk_size", 1020)
IN_DIR = config["input_dir"]

BLAST_PARAMS = config.get("blast_params",
    {
    "dust": "no",
    "xdrop_gap": "150",
    }
)
FORMATTED_BLAST_PARAMS = format_blast_params(BLAST_PARAMS)

try:
    IN_LIST = config["input_list"]
    query_genomes = []
    reference_genomes = []    

    EXPAND_FUNCTION = zip

    with open(IN_LIST) as f:
        for line in f:
            line = line.strip()
            
            if not line or line[0] == "#":
                continue
            
            genome_pair = line.split(",")
            query_genomes.append(genome_pair[0].strip())
            reference_genomes.append(genome_pair[1].strip())
        
    print("Running in paired mode")
    print("Queries: ", query_genomes)
    print("References: ", reference_genomes)

except KeyError:
    
    EXPAND_FUNCTION = itertools.product

    genomes, = glob_wildcards(os.path.join(IN_DIR, "{genome}.fna"))
    query_genomes = reference_genomes = genomes

    print("Running in all vs. all mode.\nGenomes:")
    print(genomes)


rule target:
    input: os.path.join(OUT_DIR, "anib.csv")


rule make_db:
    input: os.path.join(IN_DIR, "{reference_genome}.fna")
    output: os.path.join(OUT_DIR, "blastdb/{reference_genome}.nsq")
    params: 
        db_prefix = os.path.join(OUT_DIR, "blastdb/{reference_genome}")
    shell:
        "makeblastdb -dbtype nucl -in {input} -out {params.db_prefix}" 


rule split_query_genome:
    input: os.path.join(IN_DIR, "{query_genome}.fna")
    output: os.path.join(OUT_DIR, "split/{query_genome}.fasta")
    params:
        chunk_size = CHUNK_SIZE
    script: "scripts/split_fasta.py"


rule blast_search:
    input: query=os.path.join(OUT_DIR, "split/{query_genome}.fasta"),
           reference=os.path.join(OUT_DIR, "blastdb/{reference_genome}.nsq")
    output: os.path.join(OUT_DIR, "blast-out/{query_genome}_vs_{reference_genome}.blast.tsv")
    params: 
        db_prefix = lambda wildcards: os.path.join(OUT_DIR, f"blastdb/{wildcards.reference_genome}"),
        blast_params = FORMATTED_BLAST_PARAMS
    threads: config.get("blast_threads", 1)
    shell:
        "blastn -query {input.query} -db {params.db_prefix} -out {output} -num_threads {threads} -outfmt '6 qseqid qstart qend qlen sseqid sstart send length gaps pident nident evalue' -task blastn {params.blast_params}"
        # "blastn -query {input.query} -db {params.db_prefix} -out {output} -num_threads {threads} -outfmt '6 qseqid qstart qend qlen sseqid sstart send length gaps pident nident evalue' -dust no -penalty -1 -reward 1 -task blastn"


rule process_blast:
    input: os.path.join(OUT_DIR, "blast-out/{query_genome}_vs_{reference_genome}.blast.tsv")
    output: os.path.join(OUT_DIR, "blast-processed/{query_genome}_vs_{reference_genome}.blast.tsv")
    params: eval_threshold = config.get("eval_filter", 1.0E-15),
            coverage_threshold = config.get("coverage_filter", 0.7),
            identity_threshold = config.get("identity_filter", 0.3)
    script: "scripts/process_blast_results.py"


rule calculate_ANI:
    input: expand(os.path.join(OUT_DIR, "blast-processed/{query_genome}_vs_{reference_genome}.blast.tsv"), EXPAND_FUNCTION, query_genome = query_genomes, reference_genome = reference_genomes)
    output: os.path.join(OUT_DIR, "anib.csv")
    script: "scripts/calculate_ani.py"

