process index_target_db {
  label 'mmseqs2'
  publishDir "${params.output}/", mode: 'copy', pattern: "${targetdb.getSimpleName()}_index.tar.gz"
  publishDir "${params.output}/", mode: 'copy', pattern: 'tmp.tar.gz'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "index_target_db.log"}


  input:
  file targetdb

  output:
  tuple path("${targetdb.getSimpleName()}_index.tar.gz"), path("tmp.tar.gz"), emit: output
  file ".command.log"

  script:
  """
  echo "-----------------   Index target_db   -----------------------"
  tar -xzf ${targetdb}
  mkdir -p "${targetdb.getSimpleName()}_index"

  (cd ${targetdb.getSimpleName()}/; mmseqs createindex ${targetdb.getSimpleName()} tmp --threads ${params.threads})

  mv ${targetdb.getSimpleName()}/${targetdb.getSimpleName()}.idx ${targetdb.getSimpleName()}/${targetdb.getSimpleName()}.idx.dbtype ${targetdb.getSimpleName()}/${targetdb.getSimpleName()}.idx.index ${targetdb.getSimpleName()}_index/
  tar czf ${targetdb.getSimpleName()}_index.tar.gz ${targetdb.getSimpleName()}_index/
  mv ${targetdb.getSimpleName()}/tmp/ .
  tar czf tmp.tar.gz tmp/

  # clean-up the unzipped files
  rm -rf ${targetdb.getSimpleName()}/ ${targetdb.getSimpleName()}_index/ tmp/
  """

}
