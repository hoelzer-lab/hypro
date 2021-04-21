/* Comment section:
  - Note: 'file' needs to be replaced by 'path' in order to use emit (except when using tuple of files)
*/

process query_fasta {
  label 'query_fasta'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "query_fasta.log"}

  input:
  tuple val(name), file(prokka_out)

  output:
  path "query.fasta", emit: queryfasta
  tuple file("HyProt_loc.json"), file("HyProt_content.json"), file("gff_content.json"), emit: hyprot_dicts
  file ".command.log"

  script:
  """
  echo "----------------   Unpack prokka output   ----------------"
  tar -xzvf ${prokka_out}
  echo "----------------   Read prokka annotations   ----------------"
  query_fasta.py -gff prokka/${name}.gff -ffn prokka/${name}.ffn -m ${params.modus}
  """

}
