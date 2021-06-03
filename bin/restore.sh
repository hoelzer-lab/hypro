#!/bin/bash
{
  read
  COUNTER=0
  while IFS=$'\t\r' read ORIGINAL_ID RENAMED_ID; do
    sed -i -E "s/($RENAMED_ID)(:|_|\"|\s|$)/$ORIGINAL_ID\2/g" ${2}.*
    let COUNTER=COUNTER+1
  done
  echo "Restored contig IDs for ${COUNTER} sequences."
} < ${1}
