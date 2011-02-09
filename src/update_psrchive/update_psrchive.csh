#!/bin/tcsh -f

# Script to run on "tycho" or "newton". It will download, compile and install
# the latest version of psrchive for 32- and 64-bit computers to share.
#
# Usage: ./update_psrchive.csh [clean]
#
# If the parameter "clean" is given as the first/second argument, the build will
# begin from scratch i.e., bootstrap, configure and make clean.
#
# Author: Jonathan Khoo <Jonathan.Khoo@csiro.au>

if (`whoami` !~ "pulsar") then
	echo Please run $0 from the \"pulsar\" account
	exit
endif

set clean = 0

if ($#argv == 1) then
	if ($argv[1] =~ "clean") then
		set clean = 1
	endif
endif

set hostname = `hostname`

if ($hostname =~ "tycho") then
	set code_directory = "/pulsar/psr/cvshome/pulsar/psrchive"
else if ($hostname =~ "newton") then
	set code_directory = "/DATA/NEWTON_2/pulsar/psrchive"
else if ($hostname =~ "lagavulin") then
	set code_directory = "/psr1/cvshome/pulsar/psrchive"
else
	echo This script, $0, must be run from \"tycho\", \"newton\" or \"lagavulin\".
	exit
endif

cd $code_directory

# get the number of processors on the computer so we can hammer all of them when compiling
set num_processors = `cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'`
@ num_processors += 2

echo $hostname :: number of processors \(+ 1\): $num_processors

set pav_timestamp = `ls --full-time $PSRHOME/$LOGIN_ARCH/bin/pav | awk '{print $7}'`

if ($clean) then
	set num_steps = 7
else
	set num_steps = 4
endif

set step_count = 1

echo $hostname \($step_count/$num_steps\) :: updating code from CVS
#./update >& /dev/null
@ step_count++

if ($clean) then
	echo $hostname \($step_count/$num_steps\) :: bootstrap
	./bootstrap >& /dev/null
	@ step_count++

	echo $hostname \($step_count/$num_steps\) :: configure
	if ($hostname =~ "lagavulin") then
		./configure LOGIN_ARCH=linux PGPLOT_DIR=/psr1/packages/linux/pgplot/ PSRHOME=/psr1 --no-create --no-recursion >& /dev/null
	else
		if ($hostname =~ "newton") then
# @TODO: change either the pgplot script in psrchive or the environment defined in PSR.csh
			setenv PGPLOT_DIR /pulsar/psr/linux_64/pgplot
			setenv PGPLOT_FONT $PGPLOT_DIR/grfont.dat
		endif
		./configure >& /dev/null
	endif
	@ step_count++

	echo $hostname \($step_count/$num_steps\) :: make clean
	make clean >& /dev/null
	@ step_count++

endif

echo $hostname \($step_count/$num_steps\) :: make -j$num_processors
make -j$num_processors >& /dev/null
@ step_count++

echo $hostname \($step_count/$num_steps\) :: installing psrchive
make install >& /dev/null
@ step_count++

echo $hostname \($step_count/$num_steps\) :: checking installation

set new_pav_timestamp = `ls --full-time $PSRHOME/$LOGIN_ARCH/bin/pav | awk '{print $7}'`

if ($pav_timestamp =~ $new_pav_timestamp) then
	echo "\n\n +++ $hostname :: error with update - update manually to debug +++ \n\n"
else
	echo "\n\n --- $hostname :: automated update and installation complete --- \n\n"
endif
