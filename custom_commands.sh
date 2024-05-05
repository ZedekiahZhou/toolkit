#!/bin/bash

# get slurm job info
dirsq() {
    jname=$(grep "\-J" job* | sed -E 's/.+\-J\ //g' | tr "\n" "|" | sed 's/|$/\n/')
    if [ ${jname} ]; then sq | grep -wE "($jname)"; fi
}

dirjobID() {
    dirsq | awk '{printf $1 " "} END {print ""}'
}

myjob() {
    pkujob | awk 'BEGIN {FS="\n"; RS=""; ORS="\n\n"} /zhouzhe/'
}

myjobID() {
    myjob | grep JobId | awk -v FS="=" -v ORS=" " '{gsub(" .+", "", $2);print $2}'
}

mysq() {
    jname=$(myjob | grep JobId | sed 's/JobId=//g' | sed -E 's/\ .+//g' | tr "\n" "|" | sed 's/|$/\n/')
    if [ ${jname} ]; then sq | grep -wE "($jname)"; fi
}
