#!/bin/bash
#############################################################################################################
program_exec_string=
username=
waitdurationsecs=60
max_cns_jobs_count=

#############################################################################################################
usage()
{
cat << EOF
   USAGE for $0
	This script waits in infinite loop (with 't' seconds sleep), until 'n' number of 's' programs are running.
	If less than 'n' number of 's' programs are running, program terminates in 10 seconds.
   Options:
   -s      string that matches the programs running
   -n      maximum number of programs to be left running 
   -t      sleep time (default is 60 seconds)
   -u      username
EOF
}

while getopts "s:n:u:t:" OPTION
do
     case $OPTION in
         s)
             program_exec_string=$OPTARG
             ;;
         t)
             waitdurationsecs=$OPTARG
             ;;
         n)
             max_cns_jobs_count=$OPTARG
             ;;
         u)
             username=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done
#############################################################################################################
if [[ -z $program_exec_string ]] || [[ -z $username ]] || [[ -z $max_cns_jobs_count ]]
then
     usage
     exit 1
fi

#############################################################################################################
timestamp=`date +"%r"`
datenow=`date +"%m-%d-%y"`
job_count=`ps -u $username -f | grep $program_exec_string | grep -v grep | grep -v wait4parallel | wc -l`
echo "[$datenow $timestamp] $job_count $program_exec_string jobs running(waiting started) .."

if [ $job_count -le $max_cns_jobs_count ] 
then
	sleep 0.5
fi
 
while [ $job_count -ge $max_cns_jobs_count ]
do
	difference=`expr $job_count - $max_cns_jobs_count + 1`
	timestamp=`date +"%r"`
	datenow=`date +"%m-%d-%y"`
	echo "[$datenow $timestamp] $job_count $program_exec_string jobs already running.. waiting to finish $difference.."
	sleep $waitdurationsecs
	job_count=`ps -u bap54 -f | grep $program_exec_string | grep -v grep | grep -v wait4parallel | wc -l`
done
