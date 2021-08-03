process prokka_annotation {
  label 'prokka'
  publishDir "${params.runinfo}/${name}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "prokka_annotation.log"}

  input:
  tuple val(name), path(fasta)

  output:
  tuple val(name), path("prokka.tar.gz") , emit:output
  path ".command.log"

  script:
  """
  echo "---------------------------   Run prokka   -----------------------------"
  prokka --prefix ${name} ${params.prokka} --outdir prokka ${fasta}

  tar czf prokka.tar.gz prokka/

  # clean-up the unzipped files
  rm -rf prokka
  """

}
