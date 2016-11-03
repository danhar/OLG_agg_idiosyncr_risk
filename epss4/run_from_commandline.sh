#!/bin/bash
#
#*******************************************************************************
# Copyright (c) 2016 Daniel Harenberg - All rights reserved.
#*******************************************************************************

ulimit -f unlimited              # filesize
ulimit -d unlimited              # datasize
ulimit -s unlimited              # stacksize
ulimit -v unlimited              # virtual memory


if [[ "$HOSTNAME" == imladris ]]; then
	. /opt/intel/composerxe/bin/compilervars.sh intel64
else
	. /usr/local/intel/composerxe/bin/compilervars.sh intel64
fi

thisrun='094-testnew2'
BUILD=Build_Parallel_Optim #Optim #Debug #Parallel_Debug # Parallel_Optim #Run
projectname=epss4
export OMP_NUM_THREADS=2  # number of OpenMP threads
export OMP_STACKSIZE=16M #stack size for each OMP thread: default 4M, recommended 16M, mytest 512M
MYDIR="$( cd "$( dirname "$0" )" && pwd )"

path_to_results=~/work/Research/EPSocSec/test_results/$HOSTNAME
mkdir --parents $path_to_results/$thisrun/model_output
cd $MYDIR

cd model_input/last_results
# The next line is crazy: it finds the name of all selected calibrations (not outcommented, and without the .txt) and pipes them to tar. Note it handles also several calibs.
grep -Po '^[^!].*?(?<=/)\K[^/.]*(?=\.)' ../select_calibration_here.txt | xargs tar -cjf last_results.tar.bz2
cd ../..

tar -zcf code.tar.gz *.f90 src* run_from*
mv code.tar.gz $path_to_results/$thisrun
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod $BUILD/src_classes/*.o $BUILD/src_classes/*.mod

cp -r --parents *.f90 $BUILD src_* model_input/*.* model_input/calib* model_input/data model_input/last_results/last_results.tar.bz2 $path_to_results/$thisrun
rm model_input/last_results/last_results.tar.bz2
cd $path_to_results/$thisrun/model_input/last_results 
tar -xjf last_results.tar.bz2
rm last_results.tar.bz2
cd ../..

echo $HOSTNAME >bsub_output.txt 2>&1
echo " "  >>bsub_output.txt 2>&1

cd $BUILD
make clean >>../bsub_output.txt 2>&1
make all >>../bsub_output.txt 2>&1
cd ..

if [[ "$BUILD" == *Debug* ]]; then
    gdb --batch --eval-command=run --eval-command=bt --nw --args ./$BUILD/$projectname >>model_output/shell_output.txt 2>&1
    # module load valgrind # use either gdb or valgrind
    # valgrind --error-limit=no  ./$BUILD/$projectname >>shell_output.txt 2>&1 # not tested
else
    ./$BUILD/$projectname >model_output/shell_output.txt 2>&1
fi

if [ -d "$path_to_results/$thisrun" ]; then
    cd $path_to_results/$thisrun && rm -rf *.f90 $BUILD src_utilities src_classes
fi