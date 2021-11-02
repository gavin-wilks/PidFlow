#!/bin/bash
date

if [ $# -ne 2 ]
 then
  echo -e "\033[31m Wrong number of Arguments! \033[0m"
  exit 1
fi

energy=$1
jobid=$2
folder=19p6GeV_2019

if [ energy == 0 ]
  then 
    folder=14p5GeV_2019
elif [ energy == 1 ]
  then
    folder=19p6GeV_2019
fi

echo "${folder}"

rm -rf figures/${folder}/*
rm -rf figures/${folder}_${jobid}/*
mkdir -p figures/${folder}
mkdir -p figures/${folder}_${jobid}

root -l -b -q plotQA_RunbyRun.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_TriggerId.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Event.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Track_PID.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Track_PID_nSigma.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Track_PID_Mass2.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Track_Quality.C\(${energy},\"${jobid}\"\)
root -l -b -q plotQA_Track_Kinematics.C\(${energy},\"${jobid}\"\)

cp -av figures/${folder}/* figures/${folder}_${jobid}/.
