process create_query_db {
  label 'mmseqs2'
  publishDir "${params.output}/", mode: 'copy', pattern: 'query_db.tar.gz'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "create_query_db.log"}


  input:
  file fasta

  output:
  path "query_db.tar.gz", emit: output
  file ".command.log"

  script:
  """
  echo "-----------------   Create query_db  -----------------------"
  mmseqs createdb ${fasta} query_db
  mkdir -p tmp
  mv query_db* tmp
  mv tmp/ query_db/
  tar czf query_db.tar.gz query_db/

  # clean-up the unzipped files
  rm -rf query_db/
  """

}
