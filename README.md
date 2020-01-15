<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Language-Python3.7-green.svg)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/martinhoelzer?label=%40martinhoelzer&style=social)](https://twitter.com/martinhoelzer)

Authors: Martin H&ouml;lzer, Maximilian Arlt

# ProkkaX
Protein-coding annotation extension using additional homology searches against larger databases.

## Summary

The prokkaX tool extends common protein-coding anotations made with Prokka using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (currently available: only uniprotkb) to find homologs. For searching, mmseqs2 is utilized which currently offers a fast as well as accurate sequence comparison, to the best of our knowledge.

The tool has been tested in a conda environment (v. 4.7.11). 


## Tool Composition:

- **prokkaX.py**     main script to be called by the user
- **mmseqs2.sh**     bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting

## Requirements
ProkkaX requires the provided list of software to function properly. It is recommended to clone this repository and use a conda environment for prokkaX.

```bash
git clone https://github.com/hoelzer-lab/prokkaX.git
cd prokkaX
conda create -n prokkax python=3.7 pandas=0.25.2 mmseqs2=10.6d92c prokka=1.14.0 mygene=3.1.0
conda activate prokkax
```

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|mygene|3.1.0|Might also work for other versions|
|pandas|0.25.2|Might also work for other versions.|
|mmseqs2|10.6d92c|Install in conda environment|
|prokka (recommended)|1.14.0|used for de-novo annotation of test data of chlamydia|


## Script Usage

After installing all dependencies (see conda comand above) and cloning this repository we simply use Prokka on a test genome: 

```bash
prokka --prefix testrun --outdir run/prokka test/data/GCF_000471025.2_ASM47102v2_genomic.fna
```

followed by

```bash
scripts/prokkaX.py -i run/prokka/testrun.gff -o run/prokkax -d uniprotkb -t scripts/mmseqs2.sh -m full
```

assuming that your current working directory is the previously downloaded git repository of prokkaX. Otherwise, please adjust ``scripts/`` accordingly. 

**Arguments:**  

|Short|Long|Description|
|-----|----|-----------|
|**-h**|**--help** |show this help message and exit|
|**-i**|**--input**|Path to input gff that shall be extended|   
|**-o**|**--output**|Specify PATH to a directory. prokkaX will generate the output files to PATH.|
|**-d**|**--database**|Specifiy the target db to search for annotation extension. Current available options: [uniprotkb]|
|**-t**|**--mmseq2**|Specify the path to the mmseqs2.sh. Obligatory for execution.|
|**-m**|**--modus**|Choose the modus of prokkaX to search all hypothetical proteins (full) or leave those out which gained partial information (restricted). The dinstinction of fully un-annotated and partial annotated hypothetical proteins was observed for uniprot hits. Options: [full, restricted]|
## Output

prokkaX loads all necessary data for the extension process automatically. It stores all needed information in the ``-o/--output PATH``. In ``PATH``, it creates the following directories:

* ``db`` - stores the databases you have chosen, each in an own directory
* ``mmseqs2_output`` - to store mmseqs2 output
* ``output` - all extended files from prokka will be stored here (currently: gff, ffn, faa, gbk)
