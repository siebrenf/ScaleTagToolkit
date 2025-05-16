# ScaleBio Tagmentation Toolkit
## Version 1.1.1
* Updated conda-based workflow
  * Automated ArchR installation & usage
  * Improved conda.yml files
  * Clarified conda usage in the documentation
* Fixed fastq generation for multi-library samples
* Fixed library report generation with many samples

## Version 1.1
Initial Release as `ScaleTagToolkit`

# Old ScaleATAC releases
## Version 1.1
### samples.csv format update
The samples.csv columns `fastqIndex` and `fastqName` have been renamed to
`libIndex` and `libName` for consistency across the product. The old names are
still supported.

`libIndex` is now required when launching a run starting from a sequencer runFolder (BCLs)
with multiple libraries (PCR index). `libName` will default to
`libIndex` if no other name is specified

### Other changes

* Format updates to library structure definition `json`
* `FastqIndex` file is now specified in the library structure json
* Do not publish `tagSites.bed` (redundant with BAM)
* Tweak process resource specifications for larger runs
* Documentation updates, including a runnable example with online data

## Version 1.0
Initial Release of ScaleATAC
