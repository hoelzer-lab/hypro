/*Comment section:
  INFO: in addition to nextflow.config, the python script also contains the list of valid databases
*/

process update_prokka {
  label 'update_prokka'
  publishDir "${params.output}/${outdir}/", mode: 'copy'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.update_prokka.out', saveAs: {filename -> "update_prokka.log"}


  input:
  tuple val(name), file(prokka_output)
  tuple file(hyprotloc_dict), file(hyprotcontent_dict), file(gffcontent_dict)
  tuple file(mmseqs2_output), val(outdir)

  output:
  file "prokka_updated/*"
  path ".update_prokka.out", emit:log

  script:
  """
  tar -xzf ${prokka_output}
  mkdir -p ${prokka_output.getSimpleName()}_updated
  mkdir -p ${outdir}

  echo "----------------   Update prokka hyprots with mmseqs2   ----------------"
  update_prokka.py -ms ${mmseqs2_output} -hl ${hyprotloc_dict} -hc ${hyprotcontent_dict} -gffc ${gffcontent_dict} -i_ffn ${prokka_output.getSimpleName()}/${name}.ffn -i_faa ${prokka_output.getSimpleName()}/${name}.faa -i_gbk ${prokka_output.getSimpleName()}/${name}.gbk -m ${params.modus} -d ${params.database} -n ${name} -o prokka_updated

  mv .command.out .update_prokka.out

  # clean-up the unzipped files
  rm -rf ${prokka_output.getSimpleName()}/
  """

}
