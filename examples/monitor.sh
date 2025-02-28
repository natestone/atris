#!/usr/bin/bash
#declare current dir variable 
CWD=$(pwd)
terminal=/dev/pts/1
columns=$(stty size | awk '{ print $2}')
repl() { printf "$1"'%.s' $(eval "echo {1.."$(($2))"}"); }
#tailing last log line
#printf "#%.0s" {1..127}
date
auth=$"Sasha Banjac"
echo -e "Author: $auth"
repl = $columns
repl = $columns
echo -e "### AtRIS: Monitoring in progress..."
echo -e "### AtRIS: if you close this script, the simulation will stop"
repl = $columns
repl = $columns
echo -e 'last line of the log file:'
cd $CWD
tail -q -n 1 ./*.log 
speed=$(tail -q -n 1 ./*.log|grep -Eo '+([0-9][.][0-9]+)?'|awk '{SUM+=$1}END{print SUM}')
echo -e "### AtRIS: Net simulation speed: $speed/s"
#monitoring cpu frequency
repl = $columns
echo "cpu and ram paramters:" #please modify according to your sensors output
calc(){ awk "BEGIN { print "$*" }"; }
cpuut=$(mpstat 1 1 | awk '/^Average/ {print 100-$NF,"%"}')
memfree=$(cat /proc/meminfo | grep -e MemFree -e Buffers -e SwapFree | gawk 'BEGIN{s=0}{s+=$2}END{print s}')
memfree=$(calc $memfree/1048576)
echo -e "cpu util:     $cpuut"
freq=$(lscpu | grep MHz | head -n 1 | tr -dc '0-9')
freq=$(calc $freq/1000)
echo -e "cpu freq:     $freq MHz"
sensors | grep 'in0' | head -n 1
sensors | grep 'temp3' | head -n 1
sensors | grep 'fan1' | head -n 1
echo -e "free ram:     $memfree GB"
echo -e "free disk:    $(df -k . -h | head -n 1)"
echo -e "free disk:    $(df -k . -h | tail -n 1)"
#sensors &
repl = $columns
#start monitoring file sizes
threads=$(grep processor /proc/cpuinfo | wc -l)
echo "file sizes of most recently updated files:"
ls -lhat | head -n $threads
repl = $columns
