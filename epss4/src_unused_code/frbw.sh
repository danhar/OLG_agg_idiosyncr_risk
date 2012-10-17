#!/bin/bash
#

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

projectname="${MYDIR##*/}"

thisrun='18_50steps' # "$projectname"

mkdir --parents ~/Tempdh/frbw_output/$thisrun/model_input_old
cd $MYDIR
cp -r model_input/* ~/Tempdh/frbw_output/$thisrun/model_input_old
rm Debug_IA-32/*.o Debug_IA-32/*.mod Debug_IA-32/src_utilities/*.o Debug_IA-32/src_utilities/*.mod
mv fr-* fr-$thisrun
sed -i "s/epss\w*/$projectname/g" fr-$thisrun

ssh dharenbe@frbw.grid.uni-mannheim.de "mkdir -p $thisrun/model_output $thisrun/model_input/last_results/previous"

#scp *.f90 fr-$thisrun dharenbe@frbw.grid.uni-mannheim.de:~/$thisrun
scp -r *.f90 fr-$thisrun Debug_IA-32 Matlab src_utilities dharenbe@frbw.grid.uni-mannheim.de:~/$thisrun
scp -r model_input/*.* model_input/calib* model_input/data dharenbe@frbw.grid.uni-mannheim.de:~/$thisrun/model_input
scp model_input/last_results/*.* dharenbe@frbw.grid.uni-mannheim.de:~/$thisrun/model_input/last_results

ssh dharenbe@frbw.grid.uni-mannheim.de "cd $thisrun; qsub fr-$thisrun"