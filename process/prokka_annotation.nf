/*  Comment section:
label:        Use process labels to annotate and organize workflow processes in separate groups
              which can be referenced in the configuration file to select and configure subset
              of processes having similar computing requirements. (Or print in terminal for process tracking)
publishDir:   Publish all process output files that match the pattern into defined folder
*/

process prokka_annotation {
  label 'prokka'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "prokka_annotation.log"}

  input:
  tuple val(name), file(fasta)

  output:
  tuple val(name), path("prokka.tar.gz"), emit:output
  file ".command.log"

  script:
  """
  echo "----------------   Run prokka   ----------------"
  prokka --prefix ${name} ${params.prokka} --outdir prokka ${fasta}

  tar czf prokka.tar.gz prokka/

  # clean-up the unzipped files
  rm -rf prokka
  """

}
