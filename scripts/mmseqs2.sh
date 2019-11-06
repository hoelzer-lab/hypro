#!/bin/bash
# 22.10.19
# script serves to find hypothetical protein (hyprot) annotation from prokka
# input: query fasta with hyprots
path="$1"
QFASTA="${path}/../query.fasta"
TFASTA="$2"						#"${path}/db/uniprotkb/uniprot_sprot.fasta"
QUERYDB="${path}/query_db"
TARGETDB="$3"						#"${path}/db/uniprotkb/target_db"
RESULTDB="${path}/results_db"
TMP="${path}/tmp"
OUT="${path}/mmseq2_out.tsv"

# echo ${QUERYDB}
# echo ${TARGETDB}
# echo ${RESULTDB}
# echo ${TMP}
# echo ${OUT}
# head ${QFASTA}

mmseqs createdb ${QFASTA} ${QUERYDB}
mmseqs createdb ${TFASTA} ${TARGETDB}
mkdir -p ${TMP}
mmseqs createindex ${TARGETDB} ${TMP}
mmseqs search ${QUERYDB} ${TARGETDB} ${RESULTDB} ${TMP}
mmseqs convertalis --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits' ${QUERYDB} ${TARGETDB} ${RESULTDB} ${OUT}

head ${OUT}
echo "Generate top hit table from raw mmseq2 output..."
echo "query	target	pident	alnlen	mismatch	gapopen	qlen	qstart	qend	tstart	tend	evalue	bits" > "${path}/mmseq2_out_unique.tsv" # add header
sort -u -k1,1 ${OUT} >> "${path}/mmseq2_out_unique.tsv"
head ${path}/mmseq2_out_unique.tsv

