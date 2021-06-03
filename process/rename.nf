process rename {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "rename_ids.log"}

  input:
    tuple val(name), file(fasta)

  output:
    tuple val(name), path("${name}_renamed.fasta"), emit:renamed_contigs
    path("${name}_map.tsv"), emit:contig_map
    file ".command.log"

  script:
  """
  echo "-------------------------   Rename contig IDs   ----------------------------"
  cp ${fasta} tmp.fasta
  rename_fasta.py -i tmp.fasta -m ${name}_map.tsv -o ${name}_renamed.fasta rename
  """
}
