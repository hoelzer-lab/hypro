/*Comment section:
  INFO:   In addition to nextflow.config, the python script also contains the list of valid databases.
          If the list should be extended, think of extending it no only in the config parameter but also in update_prokka.py.
*/

process update_prokka {
  label 'update_prokka'
  publishDir "${params.output}/${name}/${outdir}/", mode: 'copy'
  publishDir "${params.runinfo}/${name}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "update_prokka.log"}


  input:
  tuple val(name), file(prokka_output), file(hyprotloc_dict), file(hyprotcontent_dict), file(gffcontent_dict), file(mmseqs2_output), val(outdir)

  output:
  file "${prokka_output.getSimpleName()}_updated/*"
  tuple val(name), path(".${name}_update_prokka.out"), emit:log
  path ".command.log"

  script:
  """
  tar -xzf ${prokka_output}
  mv prokka ${prokka_output.getSimpleName()}

  mkdir -p ${prokka_output.getSimpleName()}_updated

  echo "----------------   Update prokka hyprots with mmseqs2   ----------------"
  update_prokka.py -ms ${mmseqs2_output} -hl ${hyprotloc_dict} -hc ${hyprotcontent_dict} -gffc ${gffcontent_dict} -i_ffn ${prokka_output.getSimpleName()}/${name}.ffn -i_faa ${prokka_output.getSimpleName()}/${name}.faa -i_gbk ${prokka_output.getSimpleName()}/${name}.gbk -m ${params.modus} -d ${params.database} -n ${name} -o  ${prokka_output.getSimpleName()}_updated

  cp .command.log .${name}_update_prokka.out

  # clean-up the unzipped files
  rm -rf ${prokka_output.getSimpleName()}/
  """

}
