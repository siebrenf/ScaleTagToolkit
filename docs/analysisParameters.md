# Analysis parameters

All analysis parameters can be set in a `runParams.yml` file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can either start from an Illumina sequencer runFolder (bcl files) or a directory with (pre-generated) fastq files. Specify either
* runFolder : "path/to/runFolder" <br>
or
* fastqFolder : "path/to/fastqFolder"

where fastqFolder is a directory containing all input fastq files. See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

When starting from a sequencer run folder the workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

### Sample Information
* samples : "samples.csv"

A [file](examples/samples.csv) listing all samples in the analysis with their names, barcode sequences and optional sample settings

### Reference genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](docs/genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

### ATAC Peak Regions
By default the workflow will call ATAC peaks separately for each sample. This means that the count matrices for different samples cannot be directly merged or compared. One alternative is to use a single pre-defined set of peak regions for quantification in all samples

* peaks = "peaks.bed"

*Note* this does not apply to ArchR analysis, which does not use the workflow peak calls.

## Optional and Advanced parameters
Run `nextflow run path/to/ScaleTagToolkit --help` for a description of available options and see the example [runParams.yml](examples/runParams.yml) file.

System options (compute resource requirements, etc.) as well as all parameter defaults, are in the workflow [nextflow.config](../nextflow.config).

#### Library Structure Definition
* libStructure : "${projectDir}/references/lib-atacUniversal-a.json"

The library structure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is
* Which parts of the reads represent genomic DNA and which should be masked (e.g. the Tn5 ME sequence)

A file for our standard configuration is included in `references/lib-atacUniversal-a.json`.

