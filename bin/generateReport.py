#!/usr/bin/env python
import pandas as pd
import argparse
from os.path import exists
from pathlib import Path
from typing import Union
from reportGeneration.ReportUtil import GeneralUtils
from reportGeneration.SampleReport import SampleReport
from reportGeneration.LibraryReport import LibraryReport

def validateArgs(sampleName:Union[str,None], libName:Union[str,None], sampleSheetFn:str, resultsDir:Union[str,None], qcOutputDir:Path, libStructPath:Path) -> Union[str,None]:
    '''
    Validates:
    - @sampleSheetFn and @resultsDir exist
    - specified @sampleName, @libName exist in @sampleSheetFn if passed. 
    '''
    if (not exists(sampleSheetFn)):
        raise ValueError(f"A samplesheet does not exist at the specified path: {sampleSheetFn}")
    if (sampleName is None and libName is None):
        raise ValueError(f"Specify a -sample OR -library for which to generate a report")
    if (not exists(libStructPath)):
        raise FileNotFoundError(f"The path passed for library structure does not exist: ${libStructPath}")
    # sample argument validation path
    samplesheet = pd.read_csv(sampleSheetFn)
    availableSamples = list(samplesheet['sample'])
    if (sampleName is not None and not sampleName in availableSamples):
        raise ValueError(f"The sample provided:{sampleName} doesn't exist in specified samplesheet.\n Existing samples: {availableSamples}")
    # libName argument validation 
    availableLibs = list(samplesheet['libName'].unique())
    if (libName is not None and not libName in availableLibs):
        raise ValueError(f"The library name provided:{libName} doesn't exist in specified samplesheet\n Existing libraries: {availableLibs}")
    # results argument validation 
    if (not exists(resultsDir)):
        raise FileNotFoundError(f"The results directory provided {resultsDir} doesn't exist")
    if  str(resultsDir) != '.':
        # Validating qcDir Argument in sample report mode 
        if sampleName is not None and str(resultsDir) != '.':    
            qcDirPath = resultsDir / 'QC' / sampleName / qcOutputDir
            if not exists(qcDirPath):
                raise ValueError(f"A subdirectory named '{qcOutputDir}' within sample {sampleName}'s QC folder does not exist: {qcDirPath}")
        # Validating qcDir Argument in library report mode 
        if libName is not None:
            subset = samplesheet[samplesheet['libName'] == libName]
            associatedSamples = subset['sample']
            for sample in associatedSamples:
                qcDirPath = resultsDir / 'QC' / sample / qcOutputDir
                if not exists(qcDirPath):
                    raise ValueError(f"A subdirectory named '{qcOutputDir}' within sample {sampleName}'s QC folder does not exist: {qcDirPath}")
    return (resultsDir, qcOutputDir, libStructPath)

def resolveFastqSpecificFilePath(resultsDir: Union[str,None],libName:Union[str,None]) -> Path:
    if (str(resultsDir) == '.'):
        return Path(f"metrics.json")
    else:
        return Path(resultsDir , 'demux' , f'{libName}.demux' , 'metrics.json')

def getLibName(sampleName:str,samplesheetFn:str):
    '''
    Retrieves the 'libName' field for the 'sample' equal to @sampleName 
    from the csv at @samplesheetFn

    Assumption: values in the 'sample' column of @samplesheetFn are unique 
    '''
    sampleSheet = pd.read_csv(samplesheetFn)
    sampleData = sampleSheet.loc[sampleSheet['sample'] == sampleName].iloc[0]
    return sampleData['libName']

def createReports(sampleName:Union[str,None], samplesheetFn:str, resultsDir:str, passedLibName:Union[str,None], qcDir:Path, libStructPath:Path, makeAdvanced:bool) -> None:
    '''
    Writes Interactive HTML report(s):
    - If @sampleName is not None: for sample with @sampleName
    - If @passedLibName is not None: for library with @libName
    '''
    makeSampleReport = sampleName is not None
    makeLibReport = passedLibName is not None
    libName = passedLibName if makeLibReport else getLibName(sampleName, samplesheetFn)
    demuxMetricsPath = resolveFastqSpecificFilePath(resultsDir, libName)
    demuxJson = GeneralUtils.readJSON(demuxMetricsPath)
    if (makeSampleReport):
        reportObj = SampleReport(sampleName, resultsDir, makeAdvanced, demuxJson, qcDir, libStructPath)
        reportObj.build()
    if makeLibReport:
        reportObj = LibraryReport(libName, samplesheetFn, resultsDir, makeAdvanced, demuxJson, qcDir, libStructPath)
        reportObj.build()
    
def main():
    parser = argparse.ArgumentParser(description="An interactive HTML report generator for the ScaleBio Tag Toolkit pipeline")
    parser.add_argument("--sample", metavar="sampleName", type=str, required=False, default=None,help="The name of the sample for which a report is being generated")
    parser.add_argument("--samplesheet", metavar="sampleSheet", required=True, default=None, help="Path to the samples.csv containing information about samples and libraries")
    parser.add_argument("--results", metavar="resultsFolder", type=str, help="Path to the output folder for which the report is to be generated", required=False, default='.')
    parser.add_argument("--libName", metavar="libName", default=None, type=str,  help="libName (from passed samplesheet) for which to generate a library report")
    parser.add_argument("--advanced", default=False, action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--qcDir", default='.', required=False, help="Name of the directory in the QC output subdirectory containing the QC output to incorperate into the report")
    parser.add_argument("--libStruct", required=True, help="Library structure used to create barcode figures")
    args = parser.parse_args()
    (resultsDir, qcSubDir, libStructPath) = validateArgs(args.sample, args.libName, args.samplesheet, Path(args.results.rstrip('/')), Path(args.qcDir.rstrip('/')), Path(args.libStruct.rstrip('/')))
    createReports(args.sample, args.samplesheet, resultsDir, args.libName, qcSubDir, libStructPath, args.advanced)

if __name__ == "__main__":
    main()
