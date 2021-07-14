#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Nextflow -- Analysis Pipeline to annotate hypothetical genes from PROKKA
Authors: Kaffee-Max, Eva, Martin
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


/**************************
* INPUT CHANNELS
**************************/

// fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}



/**************************
* PROCESSES
**************************/

/* include processes that should be used outside of a sub-workflow logic */

include { rename } from './process/rename'
include { prokka_annotation } from './process/prokka_annotation'
include { restore } from './process/restore'
include { query_fasta } from './process/query_fasta'
include { download_db } from './process/download_db'
include { create_query_db ; create_query_db as create_target_db } from './process/create_query_db'
include { index_target_db } from './process/index_target_db'
include { mmseqs2 } from './process/mmseqs2'
include { update_prokka } from './process/update_prokka'
include { summary } from './process/summary'



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
      else { download_db(); db = download_db.out.db }
    }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (db_preload.exists()) { db = db_preload }
      else  { download_db(); db = download_db.out.db }
    }

  emit:
    db
}

workflow mmseqs2_dbs {
  take:
    query_db
    query_fasta

  main:
    querydb_ch = file("${params.databases_indices}/query_db.tar.gz")
    if ( !querydb_ch.exists() ) { create_query_db(query_fasta, "query_db"); querydb_ch = create_query_db.out.output }
    targetdb_ch = file("${params.databases_indices}/target_db.tar.gz")
    if ( !targetdb_ch.exists() ) { create_target_db(query_db, "target_db"); targetdb_ch = create_target_db.out.output }
    index = file("${params.databases_indices}/target_db_index.tar.gz")
    tmp = file("${params.databases_indices}/tmp.tar.gz")
    if ( !index.exists() || !tmp.exists() ) { index_target_db(targetdb_ch); targetdb_index_ch = index_target_db.out.output }
    else {targetdb_index_ch = channel.fromPath([index, tmp]).buffer( size:2 )}

  emit:
    querydb_ch
    targetdb_ch
    targetdb_index_ch
}



/**************************
* SUB-WORKFLOWS
**************************/



/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/

workflow {
      // rename contig IDs
      rename(fasta_input_ch)
      renamed_contigs = rename.out.renamed_contigs
      rename_map = rename.out.contig_map

      // run prokka annotation
      prokka_annotation(renamed_contigs)
      prokka_out_ch = prokka_annotation.out.output

      // restore original contig IDs
      restore(prokka_out_ch.join(rename_map))
      restored_prokka_contigs = restore.out.restored_contigs

      // create input fasta for mmseqs2
      query_fasta(restored_prokka_contigs)
      query_fasta_out_ch = query_fasta.out.queryfasta
      log1 = query_fasta.out.log
      hyprot_dicts_ch = query_fasta.out.hyprot_dicts

      // download query database for mmseqs2
      get_db()
      query_db = get_db.out.db

      // prepare databases for mmseqs2
      mmseqs2_dbs(query_db, query_fasta_out_ch)
      mmseqs2_querydb = mmseqs2_dbs.out.querydb_ch
      mmseqs2_targetdb = mmseqs2_dbs.out.targetdb_ch
      mmseqs2_targetdb_index = mmseqs2_dbs.out.targetdb_index_ch

      // run mmseqs2
      mmseqs2(mmseqs2_querydb, mmseqs2_targetdb, mmseqs2_targetdb_index)
      id_alninfo = mmseqs2.out.output

      // update prokka annotations
      update_prokka(restored_prokka_contigs, hyprot_dicts_ch, id_alninfo)
      log2 = update_prokka.out.log

      // produce hypro summary
      summary(log1, log2, id_alninfo)

}



/*************
* OUTPUT
*************/

workflow.onComplete {

  summary = """"""
  myFile = file("$params.output/mmseqs2_run_db${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}/hypro_summary.txt")
  myReader = myFile.newReader()
  String line
  while( line = myReader.readLine() ) {
    summary = summary + line + "\n"
  }
  myReader.close()

  log.info """
  Execution status: ${ workflow.success ? 'OK' : 'failed' }

______________________________________

\u001B[36mExecution summary\033[0m
______________________________________
$summary

Summary report:         $params.output/mmseqs2_run_db${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}/hypro_summary.txt
______________________________________

Thanks for using HYPRO!
Please cite: https://github.com/hoelzer-lab/hypro

""".stripIndent()

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
c_blue = "\u001B[36m"
c_reset = "\033[0m"
    log.info """
    ______________________________________

    ${c_blue}Workflow: HYPRO${c_reset}
    ______________________________________
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
    ______________________________________
    """.stripIndent()
}
