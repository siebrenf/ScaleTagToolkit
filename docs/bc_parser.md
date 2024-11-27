# bc_parser

`bc_parser` is a ScaleBio developed tool for combinatorial barcode handling, error-correction and demultiplexing. It supports a wide variety of library designs and sequencing configurations. It is run automatically as part of the overall workflow, but can also be used in isolation as shown here.

# Inputs
* Sequencing reads (Fastq files) `--reads` or `--read1 --read2 ...`
    * Including index read fastqs if part of the library design
    * If fastq file names follow the standard `bcl-convert` naming convention (including `_R1` for read 1, etc.), the full set of files can be given with a single option e.g. `--reads S1*.fastq.gz`).
* Library structure definition (json) `--library`
    * This defines the location and sequence-lists for all barcodes, as well as UMIs, masked sequences, etc.
    * See `references/lib-atacUniversal-a.json` for an example
* Optionally a samples.csv to demultiplex samples within a library (input fastqs) based on the first barcode-level (`--demux`)

## samples.csv
`bc_parser` uses the [samples.csv](samplesCsv.md) from the nextflow workflow input for demultiplexing multiple biological samples processed together in the ScaleBio kit. Specifically it handles splitting each library (set of fastq files) based on the `barcodes` column, corresponding to the `sample_barcode` from the library structure json (e.g. tagmentation barcode in case of the Tag Toolkit libraries)

# Command-line
A simple example run could be
```
bc_parser/bc_parser --library references/lib-atacUniversal-a.json --demux samples.csv -v --reads *fastq.gz --write-fastq --out out.demux
```

Run ``bc_parser --help` for a full list of command-line options.

# Outputs
## Reads
`bc_parser` can output the processed reads and barcodes in multiple formats. The most common case is output of transformed fastq files per sample, with error-corrected barcodes in the read-name. This is activated with `--write-fastq`

The barcodes in the read-name will be a combination of all barcoding levels (e.g. tagmentation barcode and droplet barcode) that together define a unique cell.

## Metrics
The output directory also contains an overall metrics file `metrics.json` as well as `.tsv` files with the counts of all observed barcode sequences
