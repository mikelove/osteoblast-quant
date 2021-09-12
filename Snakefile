RUNS, = glob_wildcards("fastq/{run}_R1_ALL.fastq.gz")

READS = ["R1", "R2"]

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

ANNO = "/proj/milovelab/anno"

DICT = {"129":"129S1_SvImJ", "CAST":"CAST_EiJ"}

rule all:
    input: "multiqc/multiqc_report.html"

rule salmon_index:
    input: "{ANNO}/Mus_musculus.GRCm38.v102.fa.gz"
    output: directory("{ANNO}/Mus_musculus.GRCm38.v102-salmon_1.5.2")
    shell: "{SALMON} index -p 12 -t {input} -i {output}"

rule salmon_quant_ref:
    input:
        r1 = "fastq/{sample}_R1_ALL.fastq.gz",
        r2 = "fastq/{sample}_R2_ALL.fastq.gz",
        index = "/proj/milovelab/anno/Mus_musculus.GRCm38.v102-salmon_1.5.2"
    output:
        "ref_quants/{sample}/quant.sf"
    params:
        dir = "ref_quants/{sample}"
    shell:
        "{SALMON} quant -i {input.index} -l A -p 12 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"

rule salmon_quant_diploid:
    input:
        r1 = "fastq/{strain}xB6_{sample}_R1_ALL.fastq.gz",
        r2 = "fastq/{strain}xB6_{sample}_R2_ALL.fastq.gz"
    output:
        "quants/{strain}xB6_{sample}/quant.sf"
    params:
        dir = "quants/{strain}xB6_{sample}",
        index = lambda wcs: DICT[wcs.strain]
    shell:
        "{SALMON} quant -i diploid_txomes/indices/{params.index} -l A -p 12 "
        "--numGibbsSamples 20 --thinningFactor 100 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"

rule fastqc:
    input:
        "fastq/{sample}_ALL.fastq.gz"
    output:
        "qc/{sample}/{sample}_ALL_fastqc.html"
    params:
        dir = "qc/{sample}"
    shell:
        "fastqc --quiet -t 12 --outdir {params.dir} {input}"

rule multiqc:
    input:
        expand(["ref_quants/{run}/quant.sf",
                "quants/{run}/quant.sf",
                "qc/{run}_{read}/{run}_{read}_ALL_fastqc.html"],
               run=RUNS, read=READS)
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc"
