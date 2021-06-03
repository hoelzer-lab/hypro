#!/bin/bash
# 22.10.19
# script to extend hypothetical protein (hyprot) annotation from prokka using mmseqs2
# input: query fasta with hyprots

function is_indexed(){
	#if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see gff_extend.py)
	#then
		index_files=('target_db.idx' 'target_db.idx.dbtype' 'target_db.idx.index')
		DN=$(dirname ${2})
		for file in ${index_files[@]}
		do
		# echo ${file}
		filepath=${DN}/${file}
		# echo ${filepath}
			if [ -s "${filepath}" ]
			then
				continue
				# echo "${filepath} is a file"
			else
				# echo "${filepath} corrupted"
				# echo 'Database Index is corrupted. Try to rebuild the index...'
				return 1
				echo 'Still in'
			fi
		done
		# echo ${DN}
		return 0
	# else
	# 	echo 'Unknown database type.'
	# fi

}

function created_querydb(){
	# if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see hypro.py)
	# then
		index_files=('query_db' 'query_db.dbtype' 'query_db.index' 'query_db.lookup' 'query_db_h' 'query_db_h.dbtype' 'query_db_h.index')
		DN=$(dirname ${2})
		for file in ${index_files[@]}
		do
		# echo ${file}
		filepath=${DN}/${file}
		# echo ${filepath}
			if [ -s "${filepath}" ]
			then
				continue
				# echo "${filepath} is a file"
			else
				# echo "${filepath} corrupted"
				# echo 'Database Index is corrupted. Try to rebuild the index...'
				return 1
			fi
		done
		# echo ${DN}
		return 0
	# else
	# 	echo 'Unknown database type.'
	# fi

}

function created_targetdb(){
	# if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see gff_extend.py)
	# then
		index_files=('target_db.index' 'target_db.lookup' 'target_db_h' 'target_db_h.dbtype' 'target_db_h.index')
		DN=$(dirname ${2})
		for file in ${index_files[@]}
		do
		# echo ${file}
		filepath=${DN}/${file}
		# echo ${filepath}
			if [ -s "${filepath}" ]
			then
				continue
				# echo "${filepath} is a file"
			else
				# echo "${filepath} corrupted"
				# echo 'Database Index is corrupted. Try to rebuild the index...'
				return 1
			fi
		done
		echo ${DN}
		return 0
	# else
	# 	echo 'Unknown database type.'
	# fi

}



echo "read input"
QFASTA="$1"
TFASTA="$2"						#"${path}/db/uniprotkb/uniprot_sprot.fasta"
QUERYDB="query_db"
TARGETDB="target_db"
PREFIX="results_db"
TMP="tmp"
DBTYPE="$3"
EVAL="$4"						# minimal E-value - parameter "-e"
ALEN="$5"						# minimal alignment length "--min-aln-len"
PIDENT="$6"						# minimal percent identity - "--min-seq-id"
THREADS="$7"					# number of threads to use
# SENS="$10"						# adjust sensitivity of mmseqs
RESULTDB="${RPREFIX}/db_${DBTYPE}_e${EVAL}_a${ALEN}_p${PIDENT}" # generate an individual resultsdb for every parameter setting, since the results db is not overwritten by mmseqs search
OUT="final_outs/mmseqs2_out_db_${DBTYPE}_e${EVAL}_a${ALEN}_p${PIDENT}.tsv"


echo "create content structure"
rm -rf ${PREFIX}
mkdir -p ${TMP}
mkdir -p "final_outs"
mkdir -p ${RPREFIX}


echo "prepare dbs"

created_querydb ${DBTYPE} ${QUERYDB}	# if all query_db files exist and are non-zero , skip
if [ "$?" -eq 0 ]
then
	echo 'Found a valid querydb. Skip creating query_db...'
elif [ "$?" -eq 1 ]
then
	mmseqs createdb ${QFASTA} ${QUERYDB}
fi

created_targetdb ${DBTYPE} ${TARGETDB}	# if all target_db files exist and are non-zero , skip
if [ "$?" -eq 0 ]
then
	echo 'Found a valid targetdb. Skip creating target_db...'
elif [ "$?" -eq 1 ]
then
	mmseqs createdb ${TFASTA} ${TARGETDB}
fi

is_indexed ${DBTYPE} ${TARGETDB}      # if index files of taget_db exist and non-zero, skip
if [ "$?" -eq 0 ]
then
	echo 'Found a valid index. Indexing skipped...'
elif [ "$?" -eq 1 ]
then
	mmseqs createindex ${TARGETDB} ${TMP} --threads ${THREADS}
fi


mmseqs search ${QUERYDB} ${TARGETDB} ${RESULTDB} ${TMP} --threads ${THREADS} -e ${EVAL} --min-aln-len ${ALEN} --min-seq-id ${PIDENT}
mmseqs convertalis --threads ${THREADS} --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits' ${QUERYDB} ${TARGETDB} ${RESULTDB} ${OUT}


head ${OUT}
echo "Generate unique table with highest bit scores from raw mmseq2 output..."
echo "query	target	pident	alnlen	mismatch	gapopen	qlen	qstart	qend	tstart	tend	evalue	bits" > "final_outs/mmseqs2_out_db_${DBTYPE}_e${EVAL}_a${ALEN}_p${PIDENT}_unique.tsv" # add header
sort -u -k1,1 ${OUT} >> "final_outs/mmseqs2_out_db_${DBTYPE}_e${EVAL}_a${ALEN}_p${PIDENT}_unique.tsv"
head "final_outs/mmseqs2_out_db_${DBTYPE}_e${EVAL}_a${ALEN}_p${PIDENT}_unique.tsv"
