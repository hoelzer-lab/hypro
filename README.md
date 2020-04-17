<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Language-Python3.7-green.svg)
[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer) 

Authors: Martin H&ouml;lzer, Maximilian Arlt

# HyPro
Protein-coding annotation extension using additional homology searches against larger databases.

## Summary

The HyPro tool extends common protein-coding annotations made with [Prokka](https://github.com/tseemann/prokka) using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (currently available: only [UniProtKB](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)) to find homologs. For searching, [MMseqs2](https://github.com/soedinglab/MMseqs2) is utilized which offers a fast and accurate sequence comparison.

The tool has been tested in a conda environment (v4.7.11).

## Tool Composition:

- **hypro.py**     main script to be called by the user
- **mmseqs2.sh**     bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting

## Requirements
HyPro requires the provided list of software to function properly. 

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|mygene|3.1.0|Might also work for other versions|
|pandas|0.25.2|Might also work for other versions.|
|mmseqs2|10.6d92c|Install in conda environment|
|prokka (recommended)|1.14.6|used for _de novo_ annotation of test data of chlamydia|

It is recommended to clone this repository and use a conda environment for HyPro.


```bash
git clone https://github.com/hoelzer-lab/hypro.git
cd hypro
```
Create a conda environment:

```
conda create -n hypro python=3.7 pandas=0.25.2 mmseqs2=10.6d92c prokka=1.14.6 mygene=3.1.0
conda activate hypro
```

Or simply install HyPro from the anaconda repository:

```
conda create -n hypro python=3.7
conda activate hypro
conda install hypro
```

## Script Usage

After installing all dependencies (see commands above) and cloning this repository we simply use Prokka on a test genome: 

```bash
prokka --prefix testrun --outdir run/prokka test/data/GCF_000471025.2_ASM47102v2_genomic.fna
```

and then execute the scripts from this repository:

```bash
scripts/hypro.py -i run/prokka/testrun.gff -o run/hypro -d uniprotkb -t scripts/mmseqs2.sh -m full
```
assuming that your current working directory is the previously downloaded git repository of HyPro. Otherwise, please adjust ``scripts/`` accordingly.
or in case of conda package simply do:
```bash
hypro.py -i run/prokka/testrun.gff -o run/hypro -d uniprotkb -t pathotoconda/envs/hypro/bin/mmseqs2.sh -m full
```
The mmseqs2.sh can be found in the bin directory of the conda environment where HyPro is installed.

When HyPro was already run once and the database was downloaded you can simply use the database again for other input/output data by specifiying the path to the database via ``-c``. When running the script again on the same output folder, an available database will be used again.
```bash
scripts/hypro.py -i run/prokka/testrun.gff -o run/hypro_re-use_db -t scripts/mmseqs2.sh -m full -c run/hypro/db/uniprotkb
```

**Arguments:**  

|Short|Long|Description|
|-----|----|-----------|
|**-h**|**--help** |show this help message and exit|
|**-i**|**--input**|Path to input gff that shall be extended|   
|**-o**|**--output**|Specify PATH to a directory. HyPro will generate the output files to PATH.|
|**-d**|**--database**|Specifiy the target db to search for annotation extension. Current available options: [uniprotkb]|
|**-t**|**--mmseq2**|Specify the path to the mmseqs2.sh. If using conda, the script was installed to /conda_dir/envs/my_env_name/bin/ .|
|**-m**|**--modus**|Choose the modus of HyPro to search all hypothetical proteins (full) or leave those out which gained partial annotation (restricted).|
|**-c**|**--custom-db**|Specifiy a path. HyPro will look for a db of the type defined with -d. If no database is found, HyPro will build it in the path.
partial information (restricted). The dinstinction of fully un-annotated and partial annotated hypothetical proteins was observed for uniprot hits. Options: [full, restricted]|

## Output

HyPro loads all necessary data for the extension process automatically. It stores all needed information in the ``-o/--output PATH``. In ``PATH``, it creates the following directories:

* ``db`` - stores the databases you have chosen, each in an own directory
* ``mmseqs2_output`` - to store mmseqs2 output
* ``output`` - all extended files from prokka will be stored here (currently: gff, ffn, faa, gbk)
