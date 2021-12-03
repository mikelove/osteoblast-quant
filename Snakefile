configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

ANNO = "/proj/milovelab/anno"

DICT = {"129":"129S1_SvImJ", "CAST":"CAST_EiJ"}

RUNS, = glob_wildcards("orig_fastq/{run}_R1_ALL.fastq.gz")

rule all:
     input: 
#        qc = "multiqc/multiqc_report.html"
#        ref_quants = expand("ref_quants/{run}/quant.sf", run=RUNS),
#        quants = expand("quants/{cross}_{time}/quant.sf", cross=config["crosses"], time=config["times"])
        hisat = expand("ht2_align/{cross}_{time}.bam", cross=config["crosses"], time=config["times"])

rule salmon_index_ref:
    input: "{ANNO}/Mus_musculus.GRCm38.v102.fa.gz"
    output: directory("{ANNO}/Mus_musculus.GRCm38.v102-salmon_1.5.2")
    shell: "{SALMON} index -p 12 -t {input} -i {output}"

rule salmon_quant_ref:
    input:
        r1 = "orig_fastq/{run}_R1_ALL.fastq.gz",
        r2 = "orig_fastq/{run}_R2_ALL.fastq.gz",
        index = "/proj/milovelab/anno/Mus_musculus.GRCm38.v102-salmon_1.5.2"
    output:
        "ref_quants/{run}/quant.sf"
    params:
        dir = "ref_quants/{run}"
    shell:
        "{SALMON} quant -i {input.index} -l A -p 12 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"

rule salmon_quant_diploid:
    input:
        r1 = expand("fastq/{{strain}}xB6_{{sample}}_{lane}_R1.fastq.gz", lane=config["lanes"]),
        r2 = expand("fastq/{{strain}}xB6_{{sample}}_{lane}_R2.fastq.gz", lane=config["lanes"])
    output:
        "quants/{strain}xB6_{sample}/quant.sf"
    params:
        dir = "quants/{strain}xB6_{sample}",
        index = lambda wcs: DICT[wcs.strain]
    shell:
        "{SALMON} quant -i diploid_txomes/indices/{params.index} -l A -p 12 "
        "--numBootstraps 30 --dumpEq "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"

rule fastqc:
    input:
        "orig_fastq/{run}_ALL.fastq.gz"
    output:
        "qc/{sample}/{run}_ALL_fastqc.html"
    params:
        dir = "qc/{run}"
    shell:
        "fastqc --quiet -t 12 --outdir {params.dir} {input}"

rule multiqc:
    input:
        expand(["ref_quants/{run}/quant.sf",
                "quants/{run}/quant.sf",
                "qc/{run}_{read}/{run}_{read}_ALL_fastqc.html"],
               run=RUNS, read=config["reads"])
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc"

rule hisat:
    input:
        index = "anno/grcm38_snp/genome_snp.1.ht2",
        r1 = expand("fastq/{{sample}}_{lane}_R1.fastq.gz", lane=config["lanes"]),
        r2 = expand("fastq/{{sample}}_{lane}_R2.fastq.gz", lane=config["lanes"])
    output:
        "ht2_align/{sample}.bam"
    params:
        xdir = "anno/grcm38_snp/genome_snp",
        temp_sam = "ht2_align/{sample}.unfilt.sam",
        threads = "12",
        mem = "1G",
        r1comma = lambda wildcards, input: ",".join(input.r1),
        r2comma = lambda wildcards, input: ",".join(input.r2)
    shell:
        """
        hisat2 -p {params.threads} -x {params.xdir} -1 {params.r1comma} -2 {params.r2comma} > {params.temp_sam}
        samtools view -b -q 60 {params.temp_sam} | samtools sort -@ {params.threads} -m {params.mem} -o {output}
        samtools index {output}
        rm -f {params.temp_sam}
        """
