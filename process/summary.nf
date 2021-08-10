process summary{

  input:
  tuple val(name), path(query_fasta_log)
  tuple val(name), path(update_prokka_log)
  tuple val(name), path(mmseqs2_output), val(outdir)

  output:
  path "hypro_summary.txt"

  script:
  """
  touch hypro_summary.txt
  echo "***${name}" > hypro_summary.txt
  sed -n "2p" ${query_fasta_log} >> hypro_summary.txt
  sed -n "5,6p" ${query_fasta_log} >> hypro_summary.txt
  tail -n8 ${update_prokka_log} | head -n3 >> hypro_summary.txt
  tail -n1 ${update_prokka_log} | head -n3 >> hypro_summary.txt
  """

}
