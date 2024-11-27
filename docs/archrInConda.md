## Running the preliminary ArchR analysis process using a conda environment 
#### As per the [package's own documentation](https://www.archrproject.com/index.html) the only location ArchR can be installed from is github. As an R package this means it cannot be specified as part of a Conda environment .yml. To circumvent this issue we've created a multistep process to create a conda environment containing the ArchR and its required dependencies: 
1. Ensure you have a working installation of conda installed with Python 3.8+ ([installation instructions](https://docs.conda.io/en/latest/miniconda.html#macos-installers))
2. On the command line navigate to [condaArchR](../envs/condaArchR/)
3. From within the [condaArchR](../envs/condaArchR/) directory use the following command to build a conda environment named `archrEnv` which has ArchR's dependencies installed: 
```
conda env create --quiet -f environment.yml && conda clean -a --yes
```
4. Initialize conda for your shell to ensure the newly created environment can be activated
```
conda init <shellName>
```
5. Activate your `archrEnv` environment:
```
conda activate archrEnv
```
6. With your newly created environment activated, install ArchR by running `installArchR.R`
```
Rscript installArchR.R
```
Your `archrEnv` now contains ArchR! The path to this env now needs to be specified in your params file prior to execution as `archrCondaEnv`. This path can be retrieved using the following commands. 
```
conda info --envs | grep archrEnv
```
Once this path is copied into your params file, also ensure the `runArchR` param is set to true
 
