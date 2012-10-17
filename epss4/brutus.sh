#!/bin/bash
#

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

projectname="${MYDIR##*/}"

thisrun='01b-gelsy-i13c' # "$projectname"

mkdir --parents ~/Tempdh/Brutus/$thisrun/model_input_old
cd $MYDIR
cp -r model_input/* ~/Tempdh/Brutus/$thisrun/model_input_old
rm Build_Optim/*.o Build_Optim/*.mod Build_Optim/src_utilities/*.o Build_Optim/src_utilities/*.mod
mv br-* br-$thisrun
sed -i "s/epss\w*/$projectname/g" br-$thisrun

tar -cf - *.f90 br-$thisrun Build_Optim Matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/*.* | ssh -C danielh@brutus.ethz.ch "mkdir -p $thisrun/model_output $thisrun/model_input/last_results/previous && cd $thisrun && tar -xf - ; bsub -J $thisrun < br-$thisrun"
