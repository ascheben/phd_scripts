#!/usr/bin/env bash
#Script to collect all soap logs from runJobs.sh output
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--directory)
    logpath="$2"
    shift # past argument
    ;;
    -s|--soap)
    soapname="$2"
    shift # past argument
    ;;
    -n|--name)
    logname="$2"
    shift # past argument
    ;;
esac
shift # past argument or value
done

egrep '(Query File a|Total|Paired|Single)' ${logpath}/${soapname}*  > ${logname}.mysoaplog
sample=` awk '/Query/ { print $NF }' ${logname}.mysoaplog | awk -F  "_" '{ print $1 }' `
pairs=`awk '/Pairs/  {print $(NF-1)}' ${logname}.mysoaplog`
paired=`awk '/Paired/  {print $2}' ${logname}.mysoaplog`
singled=`awk '/Singled/  {print $2}' ${logname}.mysoaplog`
elapsed=`awk '/Elapsed/ { print $NF }' ${logname}.mysoaplog`
job=`awk -F  "." '/Query/ {print $3 "." $4}' ${logname}.mysoaplog | awk -F  ":" '{print $1}'`
paste <(echo "$sample") <(echo "$pairs") <(echo "$paired") <(echo "$singled") <(echo "$elapsed") <(echo "$job") > ${logname}.tab
