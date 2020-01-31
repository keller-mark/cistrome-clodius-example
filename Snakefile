from os.path import join

SRC_DIR = "src"
DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")


rule all:
    input:
        join(PROCESSED_DIR, "chr1_127epigenomes_15observedStates.multires.mv5")

rule clodius:
    input:
        data=join(RAW_DIR, "{filename}.data.txt"),
        rowinfo=join(RAW_DIR, "{filename}.rowinfo.txt")
    output:
        join(PROCESSED_DIR, "{filename}.multires.mv5")
    shell:
        """
        clodius convert bedfile-to-multivec \
            {input.data} \
            --row-infos-filename {input.rowinfo} \
            --output-file {output} \
            --assembly hg19 \
            --starting-resolution 200 \
            --num-rows 127
        """

rule rowinfo:
    input:
        join(RAW_DIR, "{filename}.data.txt"),
        join(RAW_DIR, "{filename}.hierarchy_matrix.json")
    output:
        join(RAW_DIR, "{filename}.rowinfo.txt")
    script:
        join(SRC_DIR, "rowinfo.py")

rule split:
    input:
        join(RAW_DIR, "{filename}.hierarchy_tree.json")
    output:
        join(RAW_DIR, "{filename}.hierarchy_matrix.json")
    script:
        join(SRC_DIR, "split.py")

rule cluster:
    input:
        join(RAW_DIR, "{filename}.data.txt")
    output:
        join(RAW_DIR, "{filename}.hierarchy_tree.json")
    script:
        join(SRC_DIR, "cluster.py")


rule download:
    output:
        join(RAW_DIR, "{filename}.data.txt")
    params:
        file_url="https://github.com/Altius/epilogos/raw/master/data/{filename}.txt.gz"
    shell:
        """
        curl -L -o {output}.gz {params.file_url} && gunzip {output}.gz
        """