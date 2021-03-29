/* Comment section:
TODO:
  - Create query.fasta from prokka output (.ffn and .gff file)
    using load_gff(), is_Hyprot(), save_HyProt(), query_fasta(), load_fasta()
    from hypro script
    ==> These functions access global vars that are used by other functionalities of
        hypro
  - See mmseqs2.sh: is target_db (mmseqs2.sh input $3) an empty directory or empty file?
  - Collect process outputs in different publishDir
*/

process mmseqs2 {
  label 'mmseqs2'
  //publishDir

  input:
  file gff
  file dbfasta

  //output:

  script:
  """
  mmseqs2.sh ${params.output}/mmseqs_output/ ${dbfasta} ${params.database}/target_db ${params.database} ${params.evalue} ${params.minalnlen} ${params.pident} ${params.threads}
  """

}
