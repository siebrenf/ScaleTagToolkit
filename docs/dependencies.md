# Dependency Management

The ScaleBio Tag Toolkit workflow requires a number of dependencies to run. These include ScaleBio developed and internal executables, python libraries, etc. There are three alternative ways to provide these dependencies; select one of these, depending on what is easiest on your system, and follow the instructions below.

## Using Docker or Singularity
If your system supports [docker containers](https://www.docker.com/), this is the recommended way to handle all dependencies for the ScaleBio Tag Toolkit workflow. We provide pre-build docker containers and the workflow is setup to automatically use them.
This is enabled by adding `-profile docker` to the nextflow command-line.

If your system does not support *docker*, [singularity](https://sylabs.io/docs/) is an alternative that is enabled on many HPC clusters. Setting `-profile docker,singularity` (**no space**) will use the _singularity_ engine for all dependencies. 

See [Nextflow Containers](https://www.nextflow.io/docs/latest/container.html) for details and additional configuration options. One important point is that all input and output paths need to be available (_bind_) inside the containers. For _docker_, Nextflow should take care of that automatically, for *singularity* this requires user mounts to be enabled in the system-wide configuration (see the notes in the [Nextflow documentation](https://www.nextflow.io/docs/latest/container.html#singularity)).

## Using Conda
Another option is using the [Conda](https://docs.conda.io/en/latest) package manager to install most dependencies automatically. This is done by setting `-profile conda -with-conda`. In this case the following additional steps need to be complete
- [Manually Install ScaleBio Tools](scaleBioTools.md) on your system first.
- If running from a sequencer runFolder (.bcls) Illumina [BCL Convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) is required to be installed (and available on `$PATH`). You can find installation instructions that do not require sudo-rights [here](https://kb.10xgenomics.com/hc/en-us/articles/360001618231-How-to-troubleshoot-installing-bcl2fastq-or-bcl-convert).

## Manual Dependency installation
A list of all requirements (excluding ArchR dependencies) can be found in `envs/scaleAtac.conda.yml`, `envs/scalereport.conda.yml` and `envs/archr.conda.yml`. This can be used for manual installation in the user environment if required. All tools need to be available on `$PATH` or based in `/path/to/ScaleTagToolkit/bin/`
