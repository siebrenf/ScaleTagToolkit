# Fastq generation
We recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) from Illumina for fastq generation.

An example [samplesheet.csv](examples/samplesheet.csv) with typical options is included.

## Generate Index Read Fastqs
For ScaleBio Universal libraries the droplet barcode is sequenced in index read _I2_, while index read _I1_ is used for multiple libraries (e.g. droplet-instrument lanes). Using `bcl-convert` this can be achieved by declaring the I2 read as a `UMI` and generating index read fastqs; using `samplesheet.csv` settings:
```
CreateFastqForIndexReads,1
TrimUMI,0
OverrideCycles,Y50;I8;U16;U8Y69
```

## Retain Short Reads
We need to generate a fastq file for (short) barcode reads. In addition ATAC fragments can be quite short after trimming. We hence recommend setting 
```
MinimumTrimmedReadLength,16
MaskShortReads,16
```

## Adapter Trimming
Adapter trimming can be performed directly during fastq generation or (optionally) later on the input fastq files during the workflow. To trim during fastq generation (faster) use 
```
AdapterRead1,CTGTCTCTTATACACATCT
AdapterRead2,CTGTCTCTTATACACATCT
```

## Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing all fastq files for all samples in an analysis run. The filenames should follow the pattern `<Name>_..._<Read>_...fastq.gz`, where
* `Name` is the library name (`libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`, `I2`