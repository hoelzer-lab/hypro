/*Comment section:
  - Note: python script contains list of valid databases
  - Access input files of python script via relative paths in working dir and not by parameter flag?
*/

process update_prokka {
  label 'update_prokka'
  publishDir "${params.output}/", mode: 'copy', pattern: 'output/*'
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "update_prokka.log"}


  input:
  tuple val(name), file(prokka_output)
  tuple file(hyprotloc_dict), file(hyprotcontent_dict), file(gffcontent_dict)
  file mmseqs2_output

  output:
  file "output/*"
  file ".command.log"

  script:
  """
  echo "----------------   Unpack prokka and mmseqs2 output   ----------------"
  tar -xzvf ${prokka_output}
  tar -xzvf ${mmseqs2_output}
  mkdir -p output
  echo "----------------   Update prokka hyprots with mmseqs2   ----------------"
  update_prokka.py -ms "mmseqs2_output/final_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}_unique.tsv" -hl ${hyprotloc_dict} -hc ${hyprotcontent_dict} -gffc ${gffcontent_dict} -i_ffn "prokka/${name}.ffn" -i_faa "prokka/${name}.faa" -i_gbk "prokka/${name}.gbk" -m ${params.modus} -d ${params.database} -n ${name} -o output
  
  # clean-up the unzipped files
  rm -rf prokka mmseqs2_output
  """

}
