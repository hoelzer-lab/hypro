/*Comment section:
  - input: write target_db into .sh script instead of using it as input param?
  - @Martin: why split running mmseqs2 and filtering the output?
*/

process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/", mode: 'copy'

  input:
  file dbfasta
  file qfasta

  output:
  file "mmseqs2_output.tar.gz"

  script:
  """
  # run mmseqs2
  mmseqs2.sh . ${qfasta} ${dbfasta} target_db ${params.database} ${params.evalue} ${params.minalnlen} ${params.pident} ${params.threads}
  mkdir -p mmseqs2_output
  mv final_outs query_db* results_db target_db* mmseqs2_output
  tar czvf mmseqs2_output.tar.gz mmseqs2_output/
  """

}
