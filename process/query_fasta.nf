process query_fasta {
  label 'query_fasta'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.query_fasta.out', saveAs: {filename -> "query_fasta.log"}

  input:
  tuple val(name), path(prokka_out)

  output:
  path "query.fasta", emit: queryfasta
  tuple path("HyProt_loc.json"), path("HyProt_content.json"), path("gff_content.json"), emit: hyprot_dicts
  path ".query_fasta.out", emit: log

  script:
  """
  tar -xzf ${prokka_out}
  mv prokka ${prokka_out.getSimpleName()}

  echo "----------------   Read hypothetical proteins from prokka annotations   ----------------"
  query_fasta.py -gff ${prokka_out.getSimpleName()}/${name}.gff -ffn ${prokka_out.getSimpleName()}/${name}.ffn -m ${params.modus}

  mv .command.out .query_fasta.out

  # clean-up the unzipped files
  rm -rf ${prokka_out.getSimpleName()}
  """

}
