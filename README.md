# ScaleBio Seq Suite: Tagmentation Toolkit Workflow

This is a Nextflow workflow to run analysis of ScaleBio Tagmentation Toolkit sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (peak-count matrix, etc.), and QC reports. Optionally initial ATAC downstream analysis with ArchR can be run automatically.

## Getting started

* First install [Nextflow](http://www.nextflow.io)
* Download this workflow to your machine
* Install [dependencies](docs/dependencies.md)
* Get the bead [barcode sequence file](#barcodes)
* Launch the small pipeline [test run](#workflow-test)
* Download / configure a reference [genome](docs/genomes.md) for your samples
* Create a [samples.csv](docs/samplesCsv.md) table for your samples
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your run


## Inputs
* Sequencing reads
    * Path to the Illumina Sequencer RunFolder (bcl files)
    * If you prefer to start from fastq files, generated outside (before) this workflow, see [Fastq generation](docs/fastqGeneration.md).
* Sample Table
    * A .csv file listing all samples in the analysis with their library (PCR) index and (optional) tagmentation sample barcode sequences. See [samples.csv](docs/samplesCsv.md).
* Reference Genome
    * The workflow requires a reference genome, including a `bowtie2` index for alignment and gene annotation; see [Reference Genomes](docs/genomes.md). 
    * Pre-built reference genomes are currently available for human (grch38) and mouse (mm39). 

### Barcodes
The workflow requires information on the location and expected sequences for all cell-barcodes, see [analysisParameters.md](docs/analysisParameters.md#library-structure-definition).

*Before running any analysis* the user needs to add a file with the bead barcode sequences from the underlying droplet system software. For ATAC this file is named `737K-cratac-v1.txt.gz` and needs to be copied to `ScaleTagToolkit/references` before any analysis can be run.


## Outputs
The workflow produces alignments (`.bam` and `fragments.bed`), a cell-by-peak count-matrix (`.mtx`), QC reports and many other files. See [Outputs](docs/outputs.md) for a full list.


## Workflow Execution
### Workflow test
A small test run, with all input data stored online, can be done with 

`nextflow run /PATH/TO/ScaleTagToolkit -profile PROFILE -params-file /PATH/TO/ScaleTagToolkit/docs/examples/runParams.yml --outDir output`

See [dependencies](docs/dependencies) for the best `PROFILE` to use on your system.

### Nextflow Command-line
**Note** that `nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`.

See the [Nextflow command-line documentation](https://www.nextflow.io/docs/latest/cli.html) for the options to run `nextflow` on different systems (including HPC clusters and cloud compute).

###  Rerunning Reporting
Once output is generated, you may wish to rerun the QC filtering and HTML report generating scripts. This is particularly useful if you'd like to change the filter thresholds to include or exclude cells in the report. See [Rerunning Scripts](docs/rerunningScripts.md) for details.

## Configuration
#### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line (e.g. `--samples samples.csv`). See [analysisParameters](docs/analysisParameters.md) for details on the options.

### Quality Control and ArchR Analysis Parameters
Additional parameters used by the QC filtering process (`cellFilter`) and automated ArchR Analysis process (`archrAnalysis`) are defined in [qcAndArchR.yml](references/parameters/qcAndArchR.yml). For 
information on how best to modify see [additionalInputParams](docs/additionalInputParams.md).

### Config File
In addition to the analysis parameters, a user-specific Nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.


## Dependency Management
Different options to provide all required dependencies are described [here](docs/dependencies.md). Follow one approach there and then run nextflow with the corresponding `-profile`.

## Running in the cloud
Nextflow itself supports execution using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Nextflow tower](https://tower.nf) offers another way to manage and execute nextflow workflows online.

# Versions and Updates
See the [Change log](changelog.md)

# License
By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio.