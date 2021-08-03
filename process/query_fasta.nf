process query_fasta {
  label 'query_fasta'
  publishDir "${params.runinfo}/${name}", mode: 'copy', pattern: '.command.out', saveAs: {filename -> "query_fasta.log"}

  input:
  tuple val(name), path(prokka_out)

  output:
  tuple val(name), path("query.fasta"), emit: queryfasta
  tuple val(name), path("HyProt_loc.json"), path("HyProt_content.json"), path("gff_content.json"), emit: hyprot_dicts
  tuple val(name), path(".${name}_query_fasta.out"), emit: log
  path ".command.out"

  script:
  """
  tar -xzf ${prokka_out}
  mv prokka ${prokka_out.getSimpleName()}

  echo "----------------   Read hypothetical proteins from prokka annotations   ----------------"
  query_fasta.py -gff ${prokka_out.getSimpleName()}/${name}.gff -ffn ${prokka_out.getSimpleName()}/${name}.ffn -m ${params.modus}

  cp .command.out .${name}_query_fasta.out

  # clean-up the unzipped files
  rm -rf ${prokka_out.getSimpleName()}
  """

}
