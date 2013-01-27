#!/bin/bash
#

thisrun='027_comp_curv20'
BUILD=Build_Parallel_Optim #Optimmax #Debugmax #Parallel_Debug #    
NTHREADS=8  # number of OpenMP threads 
MYDIR="$( cd "$( dirname "$0" )" && pwd )"
projectname="${MYDIR##*/}"

mkdir --parents ~/Tempdh/Brutus/$thisrun/model_input_old
cd $MYDIR
cp -r model_input/* ~/Tempdh/Brutus/$thisrun/model_input_old
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod
mv br-* br-$thisrun
sed -i "s/epss\w*/$projectname/g" br-$thisrun

tar -cf - *.f90 br-$thisrun $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/*.* | ssh -C danielh@brutus.ethz.ch "export BUILD=$BUILD && export OMP_NUM_THREADS=$NTHREADS && mkdir -p $thisrun/model_output $thisrun/model_input/last_results/previous && cd $thisrun && tar -xf - ; bsub -J $thisrun < br-$thisrun"
