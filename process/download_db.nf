/*Comment section:
  - from template:
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_44_collection/chlamydia_gallinacea_08_1274_3/dna/Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa.gz
    gunzip Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa.gz
  - add mkdir -p target_db?
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
    tuple file("${db_out}"), path("target_db")

  script:
    if( ${params.database} == 'uniprotkb' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
      gunzip uniprot_sprot.fasta.gz
      db_out=uniprot_sprot.fasta
      mkdir -p target_db
      """

    else if( ${params.database} == 'uniref100' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
      gunzip uniref100.fasta.gz
      db_out=uniref100.fasta
      mkdir -p target_db
      """

    else if( ${params.database} == 'uniref90' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
      gunzip uniref90.fasta.gz
      db_out=uniref90.fasta
      mkdir -p target_db
      """

    else if( ${params.database} == 'uniref50' )
      """
      wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
      gunzip uniref50.fasta.gz
      db_out=uniref50.fasta
      mkdir -p target_db
      """

    else if( ${params.database} == 'pdb' )
      """
      wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
      gunzip pdb_seqres.txt.gz
      db_out=pdb_seqres.txt
      mkdir -p target_db
      """

    else
      error "Invalid database type ${params.database}"


}
