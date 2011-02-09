#!/bin/tcsh -f

# Script to update and compile PSRCHIVE on the 32- and 64-bit ATNF machines
# Author: Jonathan Khoo <Jonathan.Khoo@csiro.au>

# @TODO: add a while loop here for authentication
if (`whoami` !~ "pulsar") then
    echo Please run $0 from the \"pulsar\" account.

    #echo Enter password for the \"pulsar\" account.
    #su --login -c "/u/kho018/make_it_all.csh;" pulsar

    exit
endif

set clean = 0

if ($#argv == 1) then
	if ($argv[1] == "clean") then
		set clean = 1
	endif
endif

if ($clean) then
	ssh -f pulsar@newton "/pulsar/psr/pipeline/update_psrchive.csh clean;"
	ssh -f pulsar@tycho "/pulsar/psr/pipeline/update_psrchive.csh clean;"
	ssh -f pulsar@lagavulin "/pulsar/psr/pipeline/update_psrchive.csh clean;"
else
	ssh -f pulsar@newton "/pulsar/psr/pipeline/update_psrchive.csh;"
	ssh -f pulsar@tycho "/pulsar/psr/pipeline/update_psrchive.csh;"
	ssh -f pulsar@lagavulin "/pulsar/psr/pipeline/update_psrchive.csh;"
endif

exit
