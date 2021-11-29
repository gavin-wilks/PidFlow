#!/bin/bash
date

if [ $# -ne 3 ]
 then
  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
  exit 1
fi

jobid=$1
numjobs=$2
epdmode=$3

rm missingjobs_${jobid}.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19p6GeV_2019/OutPut/PidFlow/EpdEp/ | grep Resolution_file_19p6GeV_2019_EpdCorrections_${epdmode}_${jobid}_${j}.root`" != "Resolution_file_19p6GeV_2019_EpdCorrections_${epdmode}_${jobid}_${j}.root" ]
#if [ "`ls /gpfs01/star/scratch/gwilks3/limitedetaphi_timing/outtree${style}_${reco}/ | grep _${j}.root`" != "${Style}Tree_${j}.root" ]
then
  echo "${j}" | tee -a  missingjobs_${jobid}.log
fi
done

