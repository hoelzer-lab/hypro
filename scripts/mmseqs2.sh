#!/bin/bash
# 22.10.19
# script serves to find hypothetical protein (hyprot) annotation from prokka
# input: query fasta with hyprots

function is_indexed(){
	if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see gff_extend.py)
	then
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
	else
		echo 'Unknown database type.'
	fi

}

function created_querydb(){
	if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see gff_extend.py)
	then
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
	else
		echo 'Unknown database type.'
	fi

}

function created_resultsdb(){
	if [ ${DBTYPE} == 'uniprotkb' ] # if clause to decide for the type of db that was chosen- default is uniprotkb (see gff_extend.py)
	then
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
	else
		echo 'Unknown database type.'
	fi

}




path="$1"
QFASTA="${path}/../query.fasta"
TFASTA="$2"						#"${path}/db/uniprotkb/uniprot_sprot.fasta"
QUERYDB="${path}/query_db"
TARGETDB="$3"						#"${path}/db/uniprotkb/target_db"
RESULTDB="${path}/results_db"
TMP="${path}/tmp"
OUT="${path}/mmseq2_out.tsv"
DBTYPE="$4"

# mmseqs createdb 
mkdir -p ${TMP}

created_querydb ${DBTYPE} ${QUERYDB}	# if all query_db files exist and are non-zero , skip 
if [ "$?" -eq 0 ]
then
	echo 'Found a valid querydb. Skip creating query_db...'
elif [ "$?" -eq 1 ]
then
	mmseqs createdb ${QFASTA} ${QUERYDB}
fi


created_resultsdb ${DBTYPE} ${TARGETDB}	# if all results_db files exist and are non-zero , skip 
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
	mmseqs createindex ${TARGETDB} ${TMP}
fi

mmseqs search ${QUERYDB} ${TARGETDB} ${RESULTDB} ${TMP}
mmseqs convertalis --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits' ${QUERYDB} ${TARGETDB} ${RESULTDB} ${OUT}

head ${OUT}
echo "Generate top hit table from raw mmseq2 output..."
echo "query	target	pident	alnlen	mismatch	gapopen	qlen	qstart	qend	tstart	tend	evalue	bits" > "${path}/mmseq2_out_unique.tsv" # add header
sort -u -k1,1 ${OUT} >> "${path}/mmseq2_out_unique.tsv"
head ${path}/mmseq2_out_unique.tsv



