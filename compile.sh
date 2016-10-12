#!/bin/bash
touch .myLock
cd build
. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
make
cd ..
rm .myLock
