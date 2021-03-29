/*  Comment section:
label:        Use process labels to annotate and organize workflow processes in separate groups
              which can be referenced in the configuration file to select and configure subset
              of processes having similar computing requirements. (Or print in terminal for process tracking)
publishDir:   Publish all process output files that match the pattern into defined folder
*/

process prokka_annotation {
  label 'prokka'
  publishDir "${params.output}/prokka/", mode: 'copy'

  input:
  tuple val(name), file(fasta)

  output:
  file "${name}*"

  script:
  """
  prokka --prefix ${name} ${fasta}
  """

}
