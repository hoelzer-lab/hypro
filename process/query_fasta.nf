/* Comment section:
  - For users with no python installed: conda env in order to use python??
  - Note: 'file' needs to be replaced by 'path' to use emit (except when using tuple of files)
*/

process query_fasta {
  label 'query_fasta'
  publishDir "${params.output}/", mode: 'copy', pattern: "mmseqs2_output.tar.gz"

  input:
  tuple val(name), file(prokka_out)

  output:
  path "query.fasta", emit: queryfasta
  tuple file("HyProt_loc.json"), file("HyProt_content.json"), file("gff_content.json"), emit: hyprot_dicts

  script:
  """
  tar -xzvf ${prokka_out}
  query_fasta.py -gff prokka/${name}.gff -ffn prokka/${name}.ffn -m ${params.modus}
  """

}
