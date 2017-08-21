#!/bin/bash
#
#*******************************************************************************
# Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
#*******************************************************************************

thisrun='102-ep2ext_norpfmed2'
BUILD=Build_Parallel_Optim #Optimmax #Debugmax #Parallel_Debug # Parallel_Optim
NTHREADS=24  # number of OpenMP threads
# OMP_STACKSIZE=16M #stack size for each OMP thread: default 4M, recommended 16M, mytest 512M

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
# projectname="${MYDIR##*/}" # current directory, but careful with symlinks
# The following is robust but may take very long, so one has to wait for the rest to be executed.
# MY_IP="$(wget -q -O - checkip.dyndns.org|sed -e 's/.*Current IP Address: //' -e 's/<.*$//')"
# The alernative following one is very fast but needs VPN connection (also within ETH) and only works on GNU Linux!
MY_IP=$(hostname -I | cut -f2 -d' ')
USERNAME="$(whoami)"

#path_to_results=~/work/Research/Eqprem_idiorisk/test_results/Euler
path_to_results="/home/$USERNAME/work/Research/EPSocSec/test_results/Euler"
mkdir --parents $path_to_results/$thisrun
cd $MYDIR

cd model_input/last_results
# The next line is crazy: it finds the name of all selected calibrations (not outcommented, and without the .txt) and pipes them to tar. Note it handles also several calibs.
grep -Po '^[^!].*?(?<=/)\K[^/.]*(?=\.)' ../select_calibration_here.txt | xargs tar -cjf last_results.tar.bz2
cd ../..

# sed -i "s/epss4\w*/$projectname/g" eulerjob # only needed if projectname changes, which it usually doesn't
tar -zcf code.tar.gz *.f90 eu* src*
mv code.tar.gz $path_to_results/$thisrun
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod $BUILD/src_classes/*.o $BUILD/src_classes/*.mod

tar -cf - *.f90 eulerjob $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/last_results.tar.bz2 | ssh -C danielh@euler.ethz.ch "export BUILD=$BUILD OMP_NUM_THREADS=$NTHREADS MY_IP=$MY_IP USERNAME=$USERNAME path_to_results=$path_to_results && mkdir -p $thisrun/model_output && cd $thisrun && tar -xf - ; cd model_input/last_results; tar -xjf last_results.tar.bz2; rm last_results.tar.bz2; cd ../..; bsub -J $thisrun < eulerjob"
rm model_input/last_results/last_results.tar.bz2
