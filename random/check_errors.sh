#!/bin/env bash

#save logs copy
cp -r logs/* checked_logs/

#go to logs folder
cd logs/

#delete any files that dont have python errors
find . -type f -print0 | xargs --null grep -Z -L "Traceback" | xargs --null rm
