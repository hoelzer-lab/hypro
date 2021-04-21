#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Nextflow -- Analysis Pipeline
Author: someone@gmail.com
*/

/**************************
* META & HELP MESSAGES
**************************/

/**************************
* Help messages, user inputs & checks
**************************/

// help message
if (params.help) { exit 0, helpMSG() }

// Log infos based on user inputs
defaultMSG()

// error codes
if (params.profile) {
  exit 1, "--profile is WRONG use -profile" }

if (workflow.profile == 'standard') {
  "NO EXECUTION PROFILE SELECTED, using [-profile local,docker]" }

if (!params.fasta) {
  exit 1, "input missing, use [--fasta]"}

if (!nextflow.version.matches('20.+')) {
  println "This workflow requires Nextflow version 20.X or greater -- You are running version $nextflow.version"
  exit 1
}

/**************************
* INPUT CHANNELS
**************************/

// fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
    .view()
}

/**************************
* PROCESSES
**************************/

/* include processes that should be used outside of a sub-workflow logic */

include { module1 } from './process/module1'
include { module2 } from './process/module2'
include { prokka_annotation } from './process/prokka_annotation'
include { query_fasta } from './process/query_fasta'
include { download_db } from './process/download_db'
include { mmseqs2 } from './process/mmseqs2'
include { update_prokka } from './process/update_prokka'


/**************************
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and hpc/cloud use via params.cloudProcess.
*/

workflow get_db {
  main:

    if (!params.customdb) { db_preload = file("${params.databases}/${params.database}/${params.database}.fasta") }
    else { db_preload = file("${params.customdb}") }
    // local storage via storeDir
    if (!params.cloudProcess) {
      if ( db_preload.exists() ) { db = db_preload}
      else { download_db(); db = download_db.out }
    }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (db_preload.exists()) { db = db_preload }
      else  { download_db(); db = download_db.out }
    }

  emit: db
}

/**************************
* SUB-WORKFLOWS
**************************/

workflow subworkflow_1 {
  take:
      fasta_input_ch
      db

  main:
    module2(module1(fasta_input_ch, db))

  emit: module2.out
}


/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/

/* Comment section:
  General Notes:
    - python scripts: args structure (order of parameters, multiple inputs per flag?)
    - create subworkflows?
*/

workflow {

      /*
      if (params.fasta) {
        subworkflow_1(fasta_input_ch, query_db)
      }
      */

      // run prokka annotation
      prokka_annotation(fasta_input_ch)
      prokka_out_ch = prokka_annotation.out.output

      // create input fasta for mmseqs2
      query_fasta(prokka_out_ch)
      query_fasta_out_ch = query_fasta.out.queryfasta
      hyprot_dicts_ch = query_fasta.out.hyprot_dicts

      // download query database for mmseqs2
      get_db()
      query_db = get_db.out.db

      // run mmseqs2
      mmseqs2(query_db, query_fasta_out_ch)
      id_alninfo = mmseqs2.out.output

      // update prokka annotations
      update_prokka(prokka_out_ch, hyprot_dicts_ch, id_alninfo)

}


/*************
* --help
*************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    Workflow: HyPro

    ${c_yellow}Usage example:${c_reset}
    nextflow run wf_template --fasta '*/*.fasta'

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}            '*.fasta'  -> One sample per file
    ${c_dim}  ..change above input to csv with:${c_reset} ${c_green}--list${c_reset} true

    ${c_yellow}Parameters:${c_reset}
    --output            Name of the result folder [default: $params.output].
    --database          Specify the target db to search for annotation extension.
                        Current available options: uniprotkb, uniref50, uniref90, uniref100, pdb.
                        Note, that searching on uniref DBs will significantly extend runtime of HyPro.
                        [default: $params.database]
    --customdb          Specify a path to an existing DB. If no DB is found, HyPro will build it.
                        Requires an according -d configuration.
    --modus             Choose the modus of HyPro to search all hypothetical proteins (full) or leave
                        those out which gained partial annotation (restricted). The dinstinction of
                        fully un-annotated and partial annotated hypothetical proteins was observed for
                        uniprot annotations. Options: [full, restricted] [default: $params.modus]
    --threads           Define the number of threads to use by mmseqs indexdb, search and convertalis.
    --evalue            Include sequence matches with < e-value threshold into the profile.
                        Requires a FLOAT >= 0.0. [default: $params.evalue]
    --minalnlen         Specify the minimum alignment length as INT in range 0 to MAX aln length.
                        [default: $params.minalnlen]
    --pident            List only matches above this sequence identity for clustering.
                        Enter a FLOAT between 0 and 1.0. [default: $params.pident]

    ${c_yellow}Options:${c_reset}
    --cores             max cores per process for local use [default: $params.cores]
    --max_cores         max cores per machine for local use [default: $params.max_cores]
    --memory            max memory for local use [default: $params.memory]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.databases]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where conda environments are cached [default: $params.cachedir]

    ${c_yellow}Profiles:${c_reset}
     -profile               local,conda
                            local,docker
                            local,singularity
                            ${c_reset}
    """.stripIndent()
}

def defaultMSG() {
    log.info """
    \u001B[1;30m______________________________________\033[0m

    \u001B[36mWorkflow: TEMPLATE\033[0m
    \u001B[1;30m______________________________________\033[0m
    Profile:                $workflow.profile
    Current User:           $workflow.userName
    Nextflow-version:       $nextflow.version
    Starting time:          $nextflow.timestamp
    Workflow hash:          $workflow.commitId

        --workdir           $params.workdir
        --databases         $params.databases
        --output            $params.output
        --cores             $params.cores
        --max_cores         $params.max_cores
        --memory            $params.memory
        --cachedir          $params.cachedir
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}
