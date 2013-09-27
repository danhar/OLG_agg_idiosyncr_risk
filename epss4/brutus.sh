#!/bin/bash
#

thisrun='080-base_ep_nx4_7d'
BUILD=Build_Parallel_Optim #Optimmax #Debugmax #Parallel_Debug #    
NTHREADS=16  # number of OpenMP threads 
MYDIR="$( cd "$( dirname "$0" )" && pwd )"
projectname="${MYDIR##*/}"

mkdir --parents ~/Tempdh/Brutus/$thisrun
cd $MYDIR

cd model_input/last_results
# The next line is crazy: it finds the name of all selected calibrations (not outcommented, and without the .txt) and pipes them to tar. Note it handles also several calibs.
grep -Po '^[^!].*?(?<=/)\K[^/.]*(?=\.)' ../select_calibration_here.txt | xargs tar -cjf last_results.tar.bz2
cd ../..

tar -zcf code.tar.gz *.f90 br* src*
mv code.tar.gz ~/Tempdh/Brutus/$thisrun
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod
mv br-* br-$thisrun
sed -i "s/epss\w*/$projectname/g" br-$thisrun

tar -cf - *.f90 br-$thisrun $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/last_results.tar.bz2 | ssh -C danielh@brutus.ethz.ch "export BUILD=$BUILD && export OMP_NUM_THREADS=$NTHREADS && mkdir -p $thisrun/model_output && cd $thisrun && tar -xf - ; cd model_input/last_results; tar -xjf last_results.tar.bz2; rm last_results.tar.bz2; cd ../..; bsub -J $thisrun < br-$thisrun"
rm model_input/last_results/last_results.tar.bz2