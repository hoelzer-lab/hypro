/*Comment section:
*/

process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/", mode: 'copy', pattern: 'mmseqs2_output.tar.gz'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "mmseqs2.log"}


  input:
  file dbfasta
  file qfasta

  output:
  path "mmseqs2_output.tar.gz", emit: output
  file ".command.log"

  script:
  """
  echo "----------------   Run mmseqs2   ----------------"
  mmseqs2.sh ${qfasta} ${dbfasta} ${params.database} ${params.evalue} ${params.minalnlen} ${params.pident} ${params.threads}
  echo "----------------   Restructure output   ----------------"
  mkdir -p mmseqs2_output
  mv final_outs query_db* results_db target_db* mmseqs2_output
  tar czvf mmseqs2_output.tar.gz mmseqs2_output/
  """

}
