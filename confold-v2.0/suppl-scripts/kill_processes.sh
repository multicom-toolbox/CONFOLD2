#!/bin/bash
# Badri Adhikari, 8/2/2016
# Kill a group of processess matching a string, at once

match_string=""
user=""

usage(){
	printf "\tPURPOSE:"
	printf "\n\t Kills all processeses matching the input string"
	printf "\n\tOPTIONS:"
	printf "\n\t -u\tUser *"
	printf "\n\t -s\tString (pattern to match) *"
	printf "\n"
}

while getopts "u:s:" OPTION
do
	case $OPTION in
		u)	user=$OPTARG ;;
		s)	match_string=$OPTARG ;;
		?)	usage
			exit ;;
	esac
done

if [[ -z $match_string ]] || [[ -z $user ]]; then
     usage
     exit 1
fi

processcount=`ps -u $user aux | grep -ie $match_string | grep -v kill_processes | grep -v grep | wc -l`
if [ $processcount -gt 0 ]
then
	echo "Killing $processcount processes running that match string '$match_string'.."
	ps -u $user aux | grep -ie $match_string | grep -v kill_processes | grep -v grep | awk '{print $2}' | xargs kill -9
else
	echo "No processes that match string '$match_string'!"
fi