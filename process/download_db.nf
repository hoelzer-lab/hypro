/*Comment section:
  - What about rusure(): check for user validation before downloading very large databases?
  - Add customdb
*/

process download_db {
  label 'download_db'
  if (params.cloudProcess) {
    publishDir "${params.databases}/${params.database}/", mode: 'copy'
  }
  else {
    storeDir "${params.databases}/${params.database}/"
  }

  output:
    file "${params.database}.fasta"

  script:
    if( "${params.database}" == 'uniprotkb' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
      gunzip uniprot_sprot.fasta.gz
      mv uniprot_sprot.fasta "${params.database}.fasta"
      """

    else if( "${params.database}" == 'uniref100' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
      gunzip uniref100.fasta.gz
      """

    else if( "${params.database}" == 'uniref90' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
      gunzip uniref90.fasta.gz
      """

    else if( "${params.database}" == 'uniref50' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
      gunzip uniref50.fasta.gz
      """

    else if( "${params.database}" == 'pdb' )
      """
      wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
      gunzip pdb_seqres.txt.gz
      mv pdb_seqres.txt "${params.database}.fasta"
      """

    else
      error "Invalid database type ${params.database}"


}
