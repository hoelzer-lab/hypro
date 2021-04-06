/*Comment section:
  - input: write target_db into .sh script instead of using it as input param?
*/

process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/mmseqs2_output/final_outs/", mode: 'copy', pattern: "*.tsv"
  publishDir "${params.databases}/${params.database}/target_db/", mode: 'copy', pattern: "target_db*"

  input:
  file dbfasta
  file qfasta

  output:
  path "final_outs/*.tsv", emit: final_out
  path "query*", emit: querydb_out
  path "target*", emit: targetdb_out
  path "results_db/*", emit: results_db_out

  script:
  """
  mmseqs2.sh . ${qfasta} ${dbfasta} target_db ${params.database} ${params.evalue} ${params.minalnlen} ${params.pident} ${params.threads}
  """

}
