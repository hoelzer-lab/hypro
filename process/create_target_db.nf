process create_target_db {
  label 'mmseqs2'
  publishDir "${params.output}/", mode: 'copy', pattern: 'target_db.tar.gz'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "create_target_db.log"}


  input:
  file fasta

  output:
  path "target_db.tar.gz", emit: output
  file ".command.log"

  script:
  """
  echo "-----------------   Create target_db   -----------------------"
  mmseqs createdb ${fasta} target_db
  mkdir -p tmp
  mv target_db* tmp
  mv tmp/ target_db/
  tar czf target_db.tar.gz target_db/

  # clean-up the unzipped files
  rm -rf target_db/
  """

}
