#!/usr/bin/env bash

filelist=`find ../log_files -maxdepth 1 -type f -not -path '*/\.*' \
| grep -e "_Log.final.out" \
| sort`

for file in $filelist
do
  sampName=`echo "${file}" \
  | sed 's/\.\.\/log_files\///g ; s/_Log.final.out//g'`

  echo "cleaning Log.final.out file for ${sampName}"

  sed 's/[aA-zZ -]*[:][ ]*$//g' ${file} \
  | awk -F"|" 'BEGIN {OFS="\t"} ;
  /^\s*$/ {next;} ;{gsub("%", "Pct", $1)} ;
  {gsub(/^[ \t]+|[ \t]+$/, "", $1)};
  {gsub(/^[ \t]+|[ \t]+$/, "", $2)}; {gsub(" ", "_", $1)};
  {gsub(",", "", $1)};
  {gsub(":", "", $1)};
  {print $1,$2}' \
  | sed 's/%//g' \
  | awk -F"\t" 'NR==1{$0="ID\t'${sampName}'"RS$0} {print $0}' > \
  ../log_files/${sampName}_STAR_Stats.tab

  echo "cleaned log file: ${sampName}_STAR_Stats.tab"

done
