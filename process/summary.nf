process summary{
  publishDir "${params.output}/${outdir}/", mode:'copy', pattern: "hypro_summary.txt"

  input:
  path query_fasta_log
  path update_prokka_log
  tuple path(mmseqs2_output), val(outdir)

  output:
  path "hypro_summary.txt"

  script:
  """
  touch hypro_summary.txt
  sed -n "2p" ${query_fasta_log} > hypro_summary.txt
  sed -n "5,6p" ${query_fasta_log} >> hypro_summary.txt
  sed -n "2,4p" ${update_prokka_log} >> hypro_summary.txt
  """

}
