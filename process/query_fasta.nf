/* Comment section:
  - Note: 'file' needs to be replaced by 'path' in order to use emit (except when using tuple of files)
*/

process query_fasta {
  label 'query_fasta'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.query_fasta.out', saveAs: {filename -> "query_fasta.log"}

  input:
  tuple val(name), file(prokka_out)

  output:
  path "query.fasta", emit: queryfasta
  tuple file("HyProt_loc.json"), file("HyProt_content.json"), file("gff_content.json"), emit: hyprot_dicts
  path ".query_fasta.out", emit: log

  script:
  """
  tar -xzf ${prokka_out}

  echo "----------------   Read hypothetical proteins from prokka annotations   ----------------"
  query_fasta.py -gff ${prokka_out.getSimpleName()}/${name}.gff -ffn ${prokka_out.getSimpleName()}/${name}.ffn -m ${params.modus}

  mv .command.out .query_fasta.out

  # clean-up the unzipped files
  rm -rf ${prokka_out.getSimpleName()}
  """

}
