#!/bin/bash
MT=`expr $SGE_TASK_ID - 1`
. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
./build/Point_$MT