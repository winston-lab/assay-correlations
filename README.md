
# assay correlations

## description

A pipeline for calculating and visualizing the correlation of different assays over a set of features. Used in [our preprint](https://www.biorxiv.org/content/early/2018/06/15/347575).

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### required files

- a FASTA file of the genome

- [bedGraph](https://genome.ucsc.edu/goldenpath/help/bedgraph.html) format coverage files of the assays to be correlated
    - "stranded" assays such as TSS-seq, NET-seq, and RNA-seq should denote the strand as part of the chromosome specification, i.e. "chrI-plus" or "chrI-minus". Files with this format are automatically generated from my other pipelines and named `*-SENSE.bedgraph` or `*-ANTISENSE.bedgraph`.

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files specifying the regions over which to correlate signal
    - the BED file will be used to generate subregions as specified, i.e. given a BED file of transcript coordinates, TSS-seq signal $\pm$ 30nt from the 5' end can be correlated with RNA-seq signal over the entire transcript

## instructions
**0**. Clone this repository.

```bash
git clone https://github.com/winston-lab/assay-correlations.git
```

**1**. Create and activate the `assay_correlations` virtual environment for the pipeline using conda. This can take a while, be patient.

```bash
# navigate into the pipeline directory
cd assay-correlations

# create the assay_correlations environment
conda env create -v -f envs/default.yaml

# activate the environment
source activate assay_correlations

# to deactivate the environment
# source deactivate
```

**2**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**3**. With the `assay_correlations` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --dry-run
```

**4**. Run the pipeline, using N cores by specifying the `--cores N` flag. 

```bash
snakemake -p --cores 8
```

