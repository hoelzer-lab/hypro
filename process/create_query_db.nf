process create_query_db {
  label 'mmseqs2'
  publishDir "${params.databases_indices}/", mode: 'copy', pattern: "*${db_type}.tar.gz"
  publishDir "${params.runinfo}/${name}", mode: 'copy', pattern: ".command.log", saveAs: {pathname -> "create_${db_type}.log"}

  input:
  tuple val(name), path(fasta)
  val db_type

  output:
  tuple val(name), path("*${db_type}.tar.gz"), emit: output
  path ".command.log"

  script:

  if ("${db_type}" == 'query_db')
    """
    echo "---------------------   Create ${db_type}  ----------------------------"
    mmseqs createdb ${fasta} ${db_type}
    mkdir -p tmp
    mv ${db_type}* tmp
    mv tmp/ ${db_type}/
    tar czf ${db_type}.tar.gz ${db_type}/
    mv ${db_type}.tar.gz ${name}_${db_type}.tar.gz

    # clean-up the unzipped paths
    rm -rf ${db_type}/
    """

  else if ("${db_type}" == 'target_db')
    """
    echo "---------------------   Create ${db_type}  ----------------------------"
    mmseqs createdb ${fasta} ${db_type}
    mkdir -p tmp
    mv ${db_type}* tmp
    mv tmp/ ${db_type}/
    tar czf ${db_type}.tar.gz ${db_type}/

    # clean-up the unzipped paths
    rm -rf ${db_type}/
    """
}
