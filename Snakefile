from os.path import join

configfile: "config.yml"

SRC_DIR = "src"
DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")

FACTORS = config['factors']



rule all:
    input:
        expand(
            join(PROCESSED_DIR, "{factor}.multires.mv5"),
            factor=FACTORS
        )

# TODO

rule clodius:
    input:
        data=join(RAW_DIR, "{factor}.{cid}.bw"),
        rowinfo=join(RAW_DIR, "{factor}.{cid}.rowinfo.txt")
    output:
        join(PROCESSED_DIR, "{factor}.multires.mv5")
    shell:
        """
        clodius convert ? # TODO
        """

# TODO

rule download:
    output:
        join(RAW_DIR, "{factor}.{cid}.bw")
    params:
        file_url=(lambda w: f"http://cistrome.org/{w.cid}") # TODO
    shell:
        """
        curl -L -o {output} {params.file_url}
        """