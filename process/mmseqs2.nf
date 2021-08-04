process mmseqs2 {
  label 'mmseqs2'
  publishDir "${params.output}/${name}/mmseqs2_run_db${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}/", mode: 'copy', pattern: 'mmseqs2_outs/*'
  publishDir "${params.runinfo}/${name}", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "mmseqs2.log"}

  input:
  tuple val(name), path(query), val(db_type), path(target), path(target_index), path(tmp)

  output:
  path "mmseqs2_outs/*"
  tuple val(name), path("mmseqs2_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}_unique.tsv"), val("mmseqs2_run_db${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}"), emit: output
  path ".command.log"

  script:
  """
  echo "----------------   Prepare mmseqs2 databases   --------------------------"
  tar -xzf ${query}
  #mv query_db qdb
  mv query_db qdb
  mv qdb/* .
  rm -rf qdb/

  tar -xzf ${target}
  mv ${target.simpleName} tdb
  mv tdb/* .
  rm -rf tdb/

  tar -xzf ${target_index}
  mv ${target_index.simpleName} tdbix
  mv tdbix/* .
  rm -rf tdbix/

  tar -xzf ${tmp}


  echo "----------------   Prepare mmseqs2 output directories   ----------------"
  mkdir -p "results_db/"
  mkdir -p "mmseqs2_outs"
  RESULTDB="results_db/db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}"
  OUT="mmseqs2_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}.tsv"


  echo "----------------   Run mmseqs2   ---------------------------------------"
  # generate an individual resultsdb for every parameter setting, since the results db is not overwritten by mmseqs search
  mmseqs search query_db ${target.simpleName} \${RESULTDB} ${tmp.simpleName} --threads ${params.threads} -e ${params.evalue} --min-aln-len ${params.minalnlen} --min-seq-id ${params.pident}
  mmseqs convertalis --threads ${params.threads} --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits' query_db ${target.simpleName} \${RESULTDB} \${OUT}

  head \${OUT}
  echo "Generate unique table with highest bit scores from raw mmseq2 output..."
  # add header
  echo "query	target	pident	alnlen	mismatch	gapopen	qlen	qstart	qend	tstart	tend	evalue	bits" > "mmseqs2_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}_unique.tsv"
  sort -u -k1,1 \${OUT} >> "mmseqs2_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}_unique.tsv"
  head "mmseqs2_outs/mmseqs2_out_db_${params.database}_e${params.evalue}_a${params.minalnlen}_p${params.pident}_unique.tsv"

  #mv results_db mmseqs2_outs

  # clean-up the unzipped files
  rm -rf ${tmp}/
  rm query_db* ${target.simpleName}*
  """

}
