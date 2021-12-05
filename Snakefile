configfile: "config.json"

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

ANNO = "/proj/milovelab/anno"

DICT = {"129":"129S1_SvImJ", "CAST":"CAST_EiJ"}

RUNS, = glob_wildcards("orig_fastq/{run}_R1_ALL.fastq.gz")

rule all:
     input: 
        # ref_quants = expand("ref_quants/{run}/quant.sf", run=RUNS),
        # quants = expand("quants/{cross}_{time}/quant.sf", cross=config["crosses"], time=config["times"]),
        # hisat = expand("ht2_align/{cross}_{time}.summary", cross=config["crosses"], time=config["times"]),
        # qc = "multiqc/multiqc_report.html",
        wasp = expand("wasp_data/{cross}_genoprob.h5", cross=config["crosses"])
    

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
                "quants/{cross}_{time}/quant.sf",
                "ht2_align/{cross}_{time}.bam",
                "qc/{run}_{read}/{run}_{read}_ALL_fastqc.html"],
                run=RUNS, read=config["reads"], 
                cross=config["crosses"], time=config["times"])
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc -f -o multiqc quants/ ht2_align/ ref_quants/ qc/"

rule hisat:
    input:
        index = "anno/grcm38_snp/genome_snp.1.ht2",
        r1 = expand("fastq/{{sample}}_{lane}_R1.fastq.gz", lane=config["lanes"]),
        r2 = expand("fastq/{{sample}}_{lane}_R2.fastq.gz", lane=config["lanes"])
    output:
        align = "ht2_align/{sample}.bam",
        summary = "ht2_align/{sample}.summary"
    params:
        xdir = "anno/grcm38_snp/genome_snp",
        temp_sam = "ht2_align/{sample}.unfilt.sam",
        threads = "12",
        mem = "1G",
        r1comma = lambda wildcards, input: ",".join(input.r1),
        r2comma = lambda wildcards, input: ",".join(input.r2)
    shell:
        """
        hisat2 -p {params.threads} -x {params.xdir} \
          --summary-file {output.summary} --new-summary \
          -1 {params.r1comma} -2 {params.r2comma} > {params.temp_sam}
        samtools view -b -q 60 {params.temp_sam} | samtools sort -@ {params.threads} -m {params.mem} -o {output.align}
        samtools index {output.align}
        rm -f {params.temp_sam}
        """

def strain_to_vcf(wildcards):
    return "diploid_txomes/snps/"+DICT[wildcards.strain]+".merged.vcf.gz"

rule split_vcfs:
    input: strain_to_vcf
    output: "wasp_data/{strain}xB6_1.vcf.gz"
    params: 
        base = "wasp_data/{strain}xB6"
    shell:
        """
        while IFS= read -r line; do
          tabix -h {input} $line > {params.base}_$line.vcf;
          bgzip {params.base}_$line.vcf;
        done < wasp_data/chroms.txt
        """

rule wasp_snp2h5:
    input: "wasp_data/{cross}_1.vcf.gz"
    output:
        genoprob = "wasp_data/{cross}_genoprob.h5",
        index = "wasp_data/{cross}_snp_index.h5",
        tab = "wasp_data/{cross}_snp_tab.h5",
        hap = "wasp_data/{cross}_haps.h5"
    params:
        base = "{cross}"
    shell:
        "snp2h5 --chrom wasp_data/chromInfo.txt --format vcf "
        "--geno_prob {output.genoprob} --snp_index {output.index} "
        "--snp_tab {output.tab} --haplotype {output.hap} "
        "wasp_data/{params.base}_*.vcf.gz "

# rule find_intersecting_snps:
#     input: 
#         bam = "ht2_align/{sample}.bam",
#         index = "data/{strain}_snp_index.h5",
#         tab = "data/{strain}_snp_tab.h5",
#         hap = "data/{strain}_haps.h5"
#     output: 
#         keep = "wasp_mapping/{sample}.keep.bam",
#         remap_bam = "wasp_mapping/{sample}.to.remap.bam",
#         remap_fq1 = "wasp_mapping/{sample}.remap.fq1.gz",
#         remap_fq2 = "wasp_mapping/{sample}.remap.fq2.gz"
#     shell:
#         """
#         {MAPPING}/find_intersecting_snps.py \
#           --is_paired_end \
#           --is_sorted \
#           --output_dir wasp_mapping \
#           --snp_tab {input.tab} \
#           --snp_index {input.index} \
#           --haplotype {input.hap} \
#           --samples data/samples \
#           {input.bam}
#         """
