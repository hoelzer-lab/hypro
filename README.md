<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Language-Python3.7-green.svg)
[![Coded by: CoffeeMA2go](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/hoelzer-lab/prokkaX/commits?author=CoffeeMA2go)
[![Supervised by](https://img.shields.io/twitter/url/https/twitter.com/martinhoelzer?label=%40martinhoelzer&style=social)](https://twitter.com/martinhoelzer)

# ProkkaX
Protein-coding annotation extension using additional homology searches against larger databases.

## Summary

The prokkaX tool extends common protein-coding anotations made with Prokka using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (currently available: only uniprotkb) to find homologs. For searching, mmseqs2 is utilized which currently offers a fast as well as accurate sequence comparison, to the best of our knowledge.

The tool has been tested in a conda environment (v. 4.7.11). 


## Tool Composition:

- **gff_extend.py** main script to be called by the user
- **mmseqs2.sh**     bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting



## Requirements
ProkkaX requires the provided list of software to function properly. It is recommended to use a conda environment for prokkaX.

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|pandas|0.25.2|Might also work for other versions.|
|mmseqs2|10.6d92c|Install in conda environment|
|prokka (recommended)|1.14.0|used for de-novo annotation of test data of chlamydia|


## Script Usage

```gff_content.py -i PATH -o PATH [-d STR] -m PATH```

**optional arguments:**  

|Short|Long|Description|
|-----|----|-----------|
|**-h**|**--help** |show this help message and exit|
|**-i**|**--input**|Path to input gff that shall be extended|   
|**-o**|**--output**|Specify PATH to a directory. prokkaX will generate the output files to PATH.|
|**-d**|**--database**|Specifiy the target db to search for annotation extension. Current available options: 'uniprotkb'|
|**-m**|**--mmseq2**|Specify the path to the mmseqs2.sh. Obligatory for extension.|

## Output

prokkaX loads all necessary data for the extension automatically. For that, it stores all needed information in the **-o/--output** PATH. Here it creates a directory called

**prokkaX** - main directory, all files needed by prokkaX stored in here

In here, it will create the following subdirectories: 

**db** - stores the databases you have chosen, each in an own directory

**mmseqs2_output** - to store mmseqs2 output

**output** - all extended files from prokka will be stored here (currently: only gff file)
