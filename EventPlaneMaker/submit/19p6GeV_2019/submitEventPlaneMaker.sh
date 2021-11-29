#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/EventPlaneMaker

##########Energy Selection##########
energy=1  # 200GeV
library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/FileList/19p6GeV_2019
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19p6GeV_2019
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/PidFlow/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

# energy=2  # 27GeV
# library=SL19b
# listPath=/star/u/sunxuhit/WorkSpace/PidFlow/FileList/27GeV_2018
# outPath=/star/data01/pwg/sunxuhit/AuAu27GeV_2018
##########Energy Selection##########

##########Mode Selection##########

outDir=EpdEp

#mode=0
#inputEpdMode=0
#epdMode=0

#mode=0
#inputEpdMode=0
#epdMode=1

mode=0
inputEpdMode=1
epdMode=2



# mode=1
# outDir=ReCenterParameter

# mode=2
# outDir=ShiftParameter

# mode=3
# outDir=ShiftParameterFull

# mode=4
# outDir=Resolution

#mode=5
#outDir=ChargedFlow
##########Mode Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/PidFlow/${outDir}
mkdir -p ${outPath}/OutPut/PidFlow/${outDir}

##########Test Production##########
star-submit-template -u ie -template testEventPlaneMaker.xml -entities mode=$mode,inputEpdMode=$inputEpdMode,epdMode=$epdMode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template EventPlaneMaker_low.xml -entities mode=$mode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitEventPlaneMakerTemp.xml -entities mode=$mode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
