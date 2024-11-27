# Reference genomes
The workflow requires a set of genome reference files and settings to run. All files and settings for a genome are defined in a [genome.json](examples/genome.json) file. When launching the workflow, a reference genome is chosen by passing the path to a specific `genome.json` in the `genome` parameter (in `params.yml` or `--genome`)

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | Required | human 
bowtie_index | Path to the bowtie2 index directory and prefix | Required | `bowtie2/GRCh38` 
gtf | Path to the gene annotation | Optional; used for TSS / promoter analysis | `filteredGTF/GRCh38_transcripts.gtf` 
| macs_genome | Name of the genome for the MACS peak-caller | optional; used to adjust MACS parameters | "hs"
speciesName | Name of the species for [OrgDb](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb) | optional; used for ArchR analysis | Homo sapiens
| archrAlias | Name of the genome for [ArchR] (https://www.archrproject.com/bookdown/getting-set-up.html) | optional; used for ArchR analysis | hg38
| annotationConvention | For non-standard ArchR genomes | optional; used for ArchR analysis | ENSEMBL

* All files (`bowtie_index`, `gtf`, ...) can be either
    - An absolute path (`/path/to/genome`)
    - A relative path starting from the location of the `genome.json` file (`genes/filtered.gtf`)
    - A AWS S3 url (s3://path/to/genome)
* The ArchR settings are only needed for downstream analysis (clustering etc.) and optional

## Pre-built genomes
Pre-build reference genomes for human and mouse are available for download:
* Human: http://scale.pub.s3.amazonaws.com/genomes/grch38.tgz
* Mouse: http://scale.pub.s3.amazonaws.com/genomes/mm39.tgz
* Human / Mouse _barnyard_:  http://scale.pub.s3.amazonaws.com/genomes/grch38_mm39.tgz

Download these to your analysis server, unpack them and then use e.g.
`--genome path/to/genomes/grch38/grch38.json`

## GTF Filtering 
The workflow includes a tool to pre-filter a GTF to a subset of features (genes / transcripts). This filtered GTF can then be used for analysis:

1. Navigate to the directory with your initial `.gtf` file.
2. Create a `gtfFilter.json` by which to filter your GTF ([example](examples/gtfFilter.json))
3. Install script dependency: `pip install pybedtools`
4. Use script to filter your gtf:

`ScaleTagToolkit/bin/filter_gtf.py -gtf <path/to/genes.gtf> -j <pathToFilterJson> -o <path/to/filtered.gtf>`

5. Replace the `gtf` field of the `genome.json` with the your newly created `filtered.gtf`
