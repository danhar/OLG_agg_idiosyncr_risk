#!/bin/bash
#

thisrun='tau002test2c'
BUILD=Build_Optimmax
MYDIR="$( cd "$( dirname "$0" )" && pwd )"
projectname="${MYDIR##*/}"

mkdir --parents ~/Tempdh/Cheops/$thisrun/model_input_old
cd $MYDIR
cp -r model_input/* ~/Tempdh/Cheops/$thisrun/model_input_old
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod
mv ch-* ch-$thisrun
sed -i "s/epss\w*/$projectname/g" ch-$thisrun

tar -cf - *.f90 ch-$thisrun $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/*.* | ssh -C dharenbe@cheops.rrz.uni-koeln.de "export BUILD=$BUILD && mkdir -p $thisrun/model_output $thisrun/model_input/last_results/previous && cd $thisrun && tar -xf - ; qsub ch-$thisrun" 