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
include { download_db }    from './process/download_db'
include { prokka_annotation } from './process/prokka_annotation'
include { query_fasta } from './process/query_fasta'


/**************************
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and hpc/cloud use via params.cloudProcess.
*/

workflow get_db {
  main:
    /* TODO: check if db already exists
    db_preload = file("${params.databases}/${params.database}/*{.fasta,.txt}")
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
    */
    download_db()
    db = download_db.out
  emit: db
}

/**************************
* SUB-WORKFLOWS
**************************/

workflow subworkflow_1 {
  take:
      fasta_input_ch
      db_ch

  main:
    module2(module1(fasta_input_ch, db_ch))

  emit: module2.out
}


/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/

/* Comment section:
*/

workflow {

      /*********** Template for subworkflows ********************
      if (params.fasta) {
        subworkflow_1(fasta_input_ch, db)
      }
      ***********************************************************/

      prokka_annotation(fasta_input_ch)
      prokka_out_ch = prokka_annotation.out

      query_fasta(prokka_out_ch)
      query_fasta_out_ch = query_fasta.out
      query_fasta_out_ch.view()

      get_db()
      db = get_db.out



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

    Workflow: Template

    ${c_yellow}Usage example:${c_reset}
    nextflow run wf_template --fasta '*/*.fasta'

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}            '*.fasta'  -> one sample per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}

    ${c_yellow}Options:${c_reset}
    --cores             max cores per process for local use [default: $params.cores]
    --max_cores         max cores per machine for local use [default: $params.max_cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Parameters:${c_reset}
    --variable1             a variable [default: $params.variable1]
    --variable2             a variable [default: $params.variable2]

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
