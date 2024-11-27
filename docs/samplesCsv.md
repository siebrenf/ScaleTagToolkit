# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their library (PCR) index and sample barcode sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. When starting from a runFolder (bcl)the `libIndex` column is also required, specifying the sequencing library (PCR index) for this sample (see below). All other columns are optional. _Column names are case sensitive!_

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | Foobar-2
libIndex | Illumina library index used with this sample | S701P or `ATCGGAC` |
libName | Name for the sequencing library / fastq files; default is to use `libIndex` | S701P
barcodes | Tn5-plate wells used for this sample | 1A-2H
expectedCells | Approximate number of single cells in this sample | 50000
subsample | (Optional) Number of input reads to use for analysis of this sample| 100000000

* `sample` and `libName` should consist only of letters, numbers, dash (-) and dot (.)
* `libIndex`: Can be one sequence, a ';'-separated list of multiple sequences or a name from `references/fastqIndex.tsv`.
* When running from pre-existing fastq file input, `libName` should match the first part of the fastq file name for this sample, e.g.: `Foo1` for `Foo1_*.fastq.gz`.
* `expectedCells` is optional. If it is left out or set to 0, the number will be estimated starting from the number of cell-barcodes with over `minUniqueReads` (100; see `--qcAndArchRParams`) reads.
* `subsample` is optional. To activate subsampling, provide a target number of reads for the sample here and set `--subsample true when running the workflow`

## Demultiplexing samples within a sequencing library
During analysis the sequencing data is first split into libraries (Fastq files) based on the PCR index. These are given in the `libName` / `libIndex` columns. If multiple samples were included in one sequencing library, these are then be demultiplexed based on the tagmentation (Tn5) barcode. In that case the same library is repeated on multiple lines, once for each sample, with the specific tn5-wells for each listed in `barcodes`. E.g.

sample | libIndex | barcodes
-- | -- | --
Foo | S701P |
Bar-a | S703P | 1A-2H
Bar-b | S703P | 3A-3H

The Tn5 wells used for each sample are given in `barcodes` as either
* An individual value (`1A`)
* A range of wells (`1A-2H`)
    * Wells are sorted first by number then by letter, i.e. `1A-1H`,`2A-H`,...
    * Note that all ranges are read in column-wise order; e.g. 1A-2C, refers to 1A-1H (all of column 1) plus 2A-2C.
* A list of values or ranges, separated by semicolon (`;`) (`1A;2A-2D`) 