
import os

IN_DIR = config["input_dir"]
OUT_DIR = config["output_dir"]
CHUNK_SIZE = config.get("chunk_size", 1020)

genomes, = glob_wildcards(os.path.join(IN_DIR, "{genome}.fna"))
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
        eval_threshold = config.get("blast_eval", "1.0E-15")
    threads: config.get("blast_threads", 4)
    shell:
        "blastn -query {input.query} -db {params.db_prefix} -out {output} -num_threads {threads} -outfmt '6 qseqid qstart qend qlen sseqid sstart send length gaps pident nident evalue' -evalue {params.eval_threshold} -dust no -xdrop_gap 150 -task blastn"


rule process_blast:
    input: os.path.join(OUT_DIR, "blast-out/{query_genome}_vs_{reference_genome}.blast.tsv")
    output: os.path.join(OUT_DIR, "blast-processed/{query_genome}_vs_{reference_genome}.blast.tsv")
    script: "scripts/process_blast_results.py"


rule calculate_ANI:
    input: expand(os.path.join(OUT_DIR, "blast-processed/{query_genome}_vs_{reference_genome}.blast.tsv"), query_genome = genomes, reference_genome = genomes)
    output: os.path.join(OUT_DIR, "anib.csv")
    script: "scripts/calculate_ani.py"

