#!/bin/bash
#

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

projectname="${MYDIR##*/}"

thisrun='tau002test' # "$projectname"

mkdir --parents ~/Tempdh/Cheops/$thisrun/model_input_old
cd $MYDIR
cp -r model_input/* ~/Tempdh/Cheops/$thisrun/model_input_old
rm Build_Optim/*.o Build_Optim/*.mod Build_Optim/src_utilities/*.o Build_Optim/src_utilities/*.mod
mv ch-* ch-$thisrun
sed -i "s/epss\w*/$projectname/g" ch-$thisrun

ssh dharenbe@cheops.rrz.uni-koeln.de "mkdir -p $thisrun/model_output $thisrun/model_input/last_results/previous"

#scp *.f90 ch-$thisrun dharenbe@cheops.rrz.uni-koeln.de:~/$thisrun
scp -r *.f90 ch-$thisrun Build_Optim Matlab src_utilities dharenbe@cheops.rrz.uni-koeln.de:~/$thisrun
scp -r model_input/*.* model_input/calib* model_input/data dharenbe@cheops.rrz.uni-koeln.de:~/$thisrun/model_input
scp model_input/last_results/*.* dharenbe@cheops.rrz.uni-koeln.de:~/$thisrun/model_input/last_results

ssh dharenbe@cheops.rrz.uni-koeln.de "cd $thisrun; qsub ch-$thisrun" 