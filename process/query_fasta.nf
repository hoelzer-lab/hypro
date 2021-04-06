/* Comment section:
<<<<<<< HEAD
  - python script:
    load_gff()
    is_Hyprot()
    save_HyProt()
    query_fasta()
    load_fasta()
  - For users with no python installed: conda env in order to use python??
  - Add output channel for hypothetical proteins (HyProt_content, HyProt_loc, gff_content)
=======
  - Create query.fasta from prokka output (.ffn and .gff file)
    using load_gff(), is_Hyprot(), save_HyProt(), query_fasta(), load_fasta()
    from hypro script
  TODO:
  - For users with no python installed: conda env in order to use python??
  - Include hypothetical proteins into output channel (HyProt_content etc. from python script)
>>>>>>> e6af400d2be5dded651a8991eb63d177f013bb90
*/

process query_fasta {
  label 'query_fasta'
  publishDir "${params.output}/", mode: 'copy'

  input:
  tuple val(name), file(prokka_out)

  output:
  file "query.fasta"

  script:
  """
  tar -xzvf ${prokka_out}
  query_fasta.py -gff prokka/${name}.gff -ffn prokka/${name}.ffn -m ${params.modus}
  """

}
