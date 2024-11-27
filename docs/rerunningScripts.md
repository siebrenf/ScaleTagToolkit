### Dependency Installation
To rerun these scripts locally, their dependencies must be installed into a conda environment which will be used for execution. Using this command from the `ScaleBio/ScaleTagToolkit/`  an ArchR env named `scalereport` will be created:

`conda env create --quiet -f envs/scalereport.conda.yml && conda clean -a --yes`

Once created this environment can be activated: <br>
`conda activate scalereport`

### Rerunning scripts
The scripts run in workflow's cellFilter and report pipeline processes are `bin/atacQcFilter.py` and `bin/generateReport.py` respectively. 
The following steps can be used to rerun either script: <br>
1. Navigate to directory containing the pipeline output folder you wish to rerun the script for (the directory specified as `outDir` in your params)
2. From this directory invoke the script as an executable by providing the absolute path to the script; specifying the help flag to learn more about the arguments to pass: <br>
`/absolutePathToScaleBioRepo/ScaleTagToolkit/bin/<atacQcFilter.py|generateReport.py> --help`
3. Using the argument descriptions provided by passing the --help flag, invoke the script specifying `-results` as the pipeline output directory (`outDir`) you wish to rerun the script for. <br>

**Pipeline output for the QC and html reporting processes are found in `<outDir>/QC` and `<ourDir>/reports` respectively. To prevent overwriting the pipeline generated output an additional argument `writeDir` can be specified when running `atacQcFilter.py`**: 
- Creating new QC output and saving it to `<outDir>/QC/test`: <br> `/absolutePathToScaleBioRepo/ScaleTagToolkit/bin/atacQcFilter.py --results <outDir> --sample <sampleName> --writeDir test`

**This newly created QC output can then be used to create an HTML report by specifying the same argument for `qcDir`. The new report will be written to a child directory of `<outDir>/reports` with the same name**:
- Creating a new HTML report based on the QC output in `<outDir>/QC/test`: <br> `/absolutePathToScaleBioRepo/ScaleTagToolkit/bin/generateReport.py --results <outDir> --sample <sampleName> --qcDir test --libStruct`


### Interpreting newly generated output: 
For information on what the generated outputs represent and how they are used [review our output documentation](outputs.md)

