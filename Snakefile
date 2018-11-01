#!/usr/bin/env python

configfile: "config.yaml"

REGIONS = config["regions"]
ASSAYS = config["assays"]
CONDITIONS = config["conditions"]

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all,

rule all:
    input:
        "config.yaml",
        expand("figures/{region}_{condition}_allassays-scatterplots.svg", region=REGIONS, condition=CONDITIONS)

rule make_stranded_annotations:
    input:
        lambda wc: REGIONS[wc.region]["annotation"]
    output:
        "annotations/{region}/{region}_stranded.bed"
    shell: """
        bash scripts/makeStrandedBed.sh {input} | LC_COLLATE=C sort -k1,1 -k2,2n  > {output}
        """

rule make_modified_annotations:
    input:
        annotation = lambda wc: REGIONS[wc.region]["annotation"] if not ASSAYS[wc.assay]["stranded"] else f"annotations/{wc.region}/{wc.region}_stranded.bed",
        fasta = config["genome"]["fasta"]
    output:
        "annotations/{region}/{region}_{assay}.bed"
    params:
        upstream = lambda wc: -(REGIONS[wc.region]["include"][wc.assay]["five_end"]),
        dnstream = lambda wc: REGIONS[wc.region]["include"][wc.assay]["three_end"]
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' {input.annotation} | bedtools slop -i stdin -g <(faidx {input.fasta} -i chromsizes) -l {params.upstream} -r {params.dnstream} -s | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule map_coverage:
    input:
        annotation = lambda wc: REGIONS[wc.region]["annotation"] if (not ASSAYS[wc.assay]["stranded"] and REGIONS[wc.region]["include"][wc.assay]["whole-annotation"]) else f"annotations/{wc.region}/{wc.region}_stranded.bed" if REGIONS[wc.region]["include"][wc.assay]["whole-annotation"] else f"annotations/{wc.region}/{wc.region}_{wc.assay}.bed",
        coverage = lambda wc: ASSAYS[wc.assay]["coverage"][wc.sample]["path"]
    output:
        "scores/{region}/{region}_{assay}_{sample}.tsv"
    shell: """
        LC_COLLATE=C sort -k1,1 -k2,2n {input.annotation} | bedtools map -a stdin -b {input.coverage} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{$7=$7*1000/($3-$2); print $0}}'| sed -e 's/-plus//g;s/-minus//g' | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule average_sample_scores:
    input:
        lambda wc: expand("scores/{region}/{region}_{assay}_{sample}.tsv", sample = [k for k,v in ASSAYS[wc.assay]["coverage"].items() if v["group"]==wc.condition], region=wc.region, assay=wc.assay)
    output:
        "scores/{region}/{region}_{assay}_{condition}_allsamples.tsv"
    shell: """
        cat {input} | LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k6,6 | bedtools groupby -g 1-4,6 -c 7 -o mean | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $6, $5, "{wildcards.assay}"}}' > {output}
        """

rule cat_assay_scores:
    input:
        lambda wc: expand("scores/{region}/{region}_{assay}_{condition}_allsamples.tsv", assay=[k for k,v in ASSAYS.items() if wc.condition in [vv["group"] for kk,vv in ASSAYS[k]["coverage"].items()]], region=wc.region, condition=wc.condition)
    output:
        "scores/{region}/{region}_{condition}_allassays.tsv"
    shell: """
        cat {input} > {output}
        """

rule plot_scatter:
    input:
        table = "scores/{region}/{region}_{condition}_allassays.tsv"
    params:
        pcount = 0.01,
        anno_label = lambda wc: REGIONS[wc.region]["label"],
        cutoff_high = lambda wc: REGIONS[wc.region]["cutoff_high"],
        cutoff_low = lambda wc: REGIONS[wc.region]["cutoff_low"],
    output:
        scatter = "figures/{region}_{condition}_allassays-scatterplots.svg"
    script:
        "scripts/assay_correlations.R"

