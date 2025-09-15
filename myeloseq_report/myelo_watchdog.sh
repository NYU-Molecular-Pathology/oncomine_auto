#!/bin/sh
# --------------------------------------------------------------
# Author: Kelsey Zhu
# Description: script to start and stop watchdog daemon process
# --------------------------------------------------------------
if [ "$1" = "start" ]; then
    python3 myelo_watchdog.py >> watchdog.log &
else
    PID=$(ps -ef |grep myelo_watchdog |awk 'NR==1{print $2}')
    if [ ! -z "$PID" ];then kill $PID
    fi
fi
