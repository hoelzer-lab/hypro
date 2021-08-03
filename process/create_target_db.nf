process create_target_db {
  label 'mmseqs2'
  publishDir "${params.databases_indices}/", mode: 'copy', pattern: "*${db_type}.tar.gz"
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {pathname -> "create_${db_type}.log"}

  input:
  path fasta
  val db_type

  output:
  tuple val(name), path("*${db_type}.tar.gz"), emit: output
  path ".command.log"

  script:
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
