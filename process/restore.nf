process restore {
  publishDir "${params.output}/", mode: 'copy', pattern: "${prokka_out.getSimpleName()}.tar.gz"
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "restore_ids.log"}

  input:
    tuple val(name), file(prokka_out)
    file(map)

  output:
    tuple val(name), path("${prokka_out.getSimpleName()}.tar.gz"), emit:restored_contigs
    file ".command.log"

  script:
  """
  tar -xzf ${prokka_out}

  echo "----------------   Restore contig IDs   ----------------"
  restore.sh ${map} ${prokka_out.getSimpleName()}/${name}

  tar czf ${prokka_out.getSimpleName()}.tar.gz ${prokka_out.getSimpleName()}/

  # clean-up the unzipped files
  rm -rf ${prokka_out.getSimpleName()}/
  """
}
