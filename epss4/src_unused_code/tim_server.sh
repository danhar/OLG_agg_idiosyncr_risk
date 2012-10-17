#!/bin/bash
#

thisrun='test'

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

projectname="${MYDIR##*/}"

mkdir ~/Tempdh/tim_output/$thisrun
cd $MYDIR # Need this coz when file is clicked from desktop, PWD is wrong

ssh tim-server.vwl.uni-mannheim.de "mkdir $thisrun; mkdir $thisrun/model_output"

scp *.f90 tim-server.vwl.uni-mannheim.de:~/$thisrun
scp -r Debug_IA-32 model_input src_parameters src_utilities tim-server.vwl.uni-mannheim.de:~/$thisrun

ssh tim-server.vwl.uni-mannheim.de "cd $thisrun/Debug_IA-32; sed -i 's/-lmkl_lapack95 /-lmkl_lapack95_lp64 /g' objects.mk; make all; cd ..; ./Debug_IA-32/$projectname &>konsole_out.txt; scp -r konsole_out.txt model_output model_input makrowi7.vwl.uni-mannheim.de:~/Tempdh/tim_output/$thisrun"
