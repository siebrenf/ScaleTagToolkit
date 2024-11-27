# Additional Input Parameter Specification 

### In addition to [runParams.yml](examples/runParams.yml) which contains the majority of configuration parameters, an additional .yml ([`qcAndArchR.yml`](../references/parameters/qcAndArchR.yml)) is used to store parameters used for Quality Control and Automated ArchR analysis.  

<br>

## Modifying QC and ArchR params 
By default [`qcAndArchR.yml`](../references/parameters/qcAndArchR.yml) is specified as the yml to be parsed to access these QC and ArchR params in the [nextflow.config](../nextflow.config). If you wish to modify these params the following steps are recommended: <br> 
1. Take a look at [`qcAndArchR.yml`](../references/parameters/qcAndArchR.yml) for a description of each parameter  
2. Make a copy of [`qcAndArchR.yml`](../references/parameters/qcAndArchR.yml) to a directory of your choice: <br>`cp ScaleBio/ScaleTagToolkit/references/parameters/qcAndArchR.yml <pathToYmlCopy>`
3. **Modify the copy at `pathToYmlCopy` to suit your needs
4. To use your newly defined .yml add: <br> `qcAndArchRParams: <pathToYmlCopy>` in your [runParams.yml](examples/runParams.yml) (the yml being passed to `-params-file` argument)
 
**Warning all parameters specified in [`qcAndArchR.yml`](../references/parameters/qcAndArchR.yml) are required and must be defined in your new .yml.


