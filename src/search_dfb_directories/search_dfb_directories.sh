#!/bin/bash

#
# Search all DFB directories ($DFB1, $DFB2, $DFB3, $DFB4, and $DFB5)
# for the target filename (and pulsar name) and prints out its full path.
#
# Usage:
#   $ search_dfb_directories <filename> <pulsar name>
#
# Author: Jonathan Khoo
# Date: 30.08.10
#
# 1. You only need pulsar name -> get rid of filename

DFB1_DIRECTORY="/pulsar/archive06/DFB"
DFB2_DIRECTORY="/pulsar/archive12/DFB"
DFB3_DIRECTORY="/pulsar/archive14/DFB"
DFB4_DIRECTORY="/pulsar/archive18/DFB"
DFB5_DIRECTORY="/pulsar/archive19/DFB"

DFB_DIRECTORIES=($DFB1_DIRECTORY $DFB2_DIRECTORY $DFB3_DIRECTORY $DFB4_DIRECTORY $DFB5_DIRECTORY)

# TODO: LOOKING HERE
FILENAME=$1
PULSAR_NAME=$2

# set a new counter (make it 0)

for directory in ${DFB_DIRECTORIES[*]}
do

  # 1. count number of files in directory 
  # 2. add that to your counter


  # result=`find ${directory}/$PULSAR_NAME -name "${FILENAME}" | wc -l | awk '{if ($1 == 0) print "0"; else print "1"}'`
  # echo ${directory}/${PULSAR_NAME}/${FILENAME}
done

# output your final counter

