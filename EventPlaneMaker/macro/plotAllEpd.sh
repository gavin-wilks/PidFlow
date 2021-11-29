#!/bin/bash
date

if [ $# -ne 3 ]
 then
  echo -e "\033[31m Wrong number of Arguments! \033[0m"
  exit 1
fi

energy=$1
epdMode=$2
jobid=$3

folder=19p6GeV_2019

if [ energy == 0 ]
  then 
    folder=14p5GeV_2019
elif [ energy == 1 ]
  then
    folder=19p6GeV_2019
fi

echo "${folder}"

mkdir -p figures/${folder}_${epdMode}_${jobid}

root -l -b -q plotEpdEp.C\(${energy},${epdMode},\"${jobid}\"\)
root -l -b -q plotEpdEpResolution.C\(${energy},${epdMode},\"${jobid}\"\)

cp -av figures/*${jobid}.pdf figures/${folder}_${epdMode}_${jobid}/.
