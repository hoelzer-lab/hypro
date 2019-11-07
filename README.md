# ProkkaX
Protein-coding annotation extension using additional homology searches against larger databases.

---------------------

## Summary

The prokkaX tool extends common protein-coding anotations made with Prokka using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (currently available: only uniprotkb) to find homologs. For searching, mmseqs2 is utilized which currently offers a fast as well as accurate sequence comparison, to the best of our knowledge.

The tool has been tested in an conda environment (v. 4.7.11). 

-----------
## Tool Composition:

- **gff_extend.py** main script to be called by the user
- **mmseqs2.sh**     bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting


----
## Requirements
ProkkaX requires the provided list of software to function properly. It is recommended to use a conda environment for prokkaX.

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|pandas|0.25.2|Might also work for ower versions.|
|mmseqs2|10.6d92c|Install in conda environment|
|prokka (recommended)|1.14.0|used for de-novo annotation of test data of chlamydia|
----

## Script Usage

```gff_content.py -i PATH -o PATH [-d STR] -m PATH```

**optional arguments:**  

|Short|Long|Description|
|-----|----|-----------|
|**[-h]**|**[--help]** |show this help message and exit|
|**-i**|**--input**|Path to input gff that shall be extended|   
|**-o**|**--output**|Specify PATH to a directory. prokkaX will generate the output files to PATH.|
|**[-d]**|**[--database]**|Specifiy the target db to search in for annotation extension. Current available options: 'uniprotkb'|
|**-m**|**--mmseq2**|Specify the path to the mmseqs2.sh. Obligatory for extension.|

## Output

prokkaX loads all necessary data for the extension automatically. For that, it stores all needed information in the **-o/--output** PATH. Here it creates a directory called

**prokkaX** - main directory, all files needed by prokkaX stored in here

In here, it will create the following subdirectories:
**db** - stores the database of choice in here
**mmseqs2_output** - to store mmseqs2 output
**output** - all extended files from prokka will be stored here (currently: only gff file)
