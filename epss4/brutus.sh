#!/bin/bash
#

thisrun='090-base_pidel_05'
BUILD=Build_Parallel_Optim #Optimmax #Debugmax #Parallel_Debug # Parallel_Optim    
NTHREADS=16  # number of OpenMP threads 
MYDIR="$( cd "$( dirname "$0" )" && pwd )"
projectname="${MYDIR##*/}"
MY_IP="$(wget -q -O - checkip.dyndns.org|sed -e 's/.*Current IP Address: //' -e 's/<.*$//')"
USERNAME="$(whoami)"

#path_to_results=~/Tempdh/Brutus
#path_to_results=~/work/Research/Eqprem_idiorisk/test_results/Brutus
path_to_results=~/work/Research/EPSocSec/test_results/Brutus
mkdir --parents $path_to_results/$thisrun
cd $MYDIR

cd model_input/last_results
# The next line is crazy: it finds the name of all selected calibrations (not outcommented, and without the .txt) and pipes them to tar. Note it handles also several calibs.
grep -Po '^[^!].*?(?<=/)\K[^/.]*(?=\.)' ../select_calibration_here.txt | xargs tar -cjf last_results.tar.bz2
cd ../..

mv br-* br-$thisrun
sed -i "s/epss4\w*/$projectname/g" br-$thisrun
tar -zcf code.tar.gz *.f90 br* src*
mv code.tar.gz $path_to_results/$thisrun
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod $BUILD/src_classes/*.o $BUILD/src_classes/*.mod

tar -cf - *.f90 br-$thisrun $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/last_results.tar.bz2 | ssh -C danielh@brutus.ethz.ch "export BUILD=$BUILD OMP_NUM_THREADS=$NTHREADS MY_IP=$MY_IP USERNAME=$USERNAME path_to_results=$path_to_results && mkdir -p $thisrun/model_output && cd $thisrun && tar -xf - ; cd model_input/last_results; tar -xjf last_results.tar.bz2; rm last_results.tar.bz2; cd ../..; bsub -J $thisrun < br-$thisrun"
rm model_input/last_results/last_results.tar.bz2