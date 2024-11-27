# ScaleBio tools

In addition to third-party and open-source software the workflow also uses several tools developed by ScaleBio:

- bc_parser
	- Extracts and error corrects cell-barcodes and UMIs from the original (input) fastq files
	- Splits (demultiplexes) the input fastq files into sample fastq files based on cell-barcodes (tagmentation wells)
	- Barcode and read-level metrics
- sc_dedup
	- BAM deduplication aware of the multi-level combinatorial cell-barcodes
	- Cell-barcode metrics
- sc_counter
	- Cell-by-peak count matrix generation, used for peak and TSS counts.
	- Peak metrics

### Installation
These tools are included in the `scaleAtac` docker container image; when running the workflow in the recommended configuration with _docker_ (or _singularity_, _podman_ etc.), they will be automatically available.

 These tool are however currently not available through Conda, so if running without docker, they need to be installed first. A download script to get static pre-compiled binaries for linux (x86_64) is provided in `envs/download-scale-tools.sh`. It will install the binaries in the `bin` directory inside the nextflow workflow, from where they will be available to the workflow.
