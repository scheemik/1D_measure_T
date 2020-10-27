#!/bin/bash
# A bash script to run the Dedalus python code
# Takes in optional arguments:
#	$ sh run.sh -n <name of experiment>         Default: current datetime
#							-c <cores>                      Default: 1
#							-r <run simulation>             Default: False
#             -m <merge h5 files>             Default: False
#							-o <post-process data>          Default: False
#							-p <plot simulation>            Default: False
#             -g <create gif>                 Default: False
#             -v <create video, mp4>          Default: False

# Current datetime
DATETIME=`date +"%Y-%m-%d_%Hh%M"`

# Having a ":" after a flag means an option is required to invoke that flag
while getopts n:c:rmopgv option
do
	case "${option}"
		in
		n) NAME=${OPTARG};;
		c) CORES=${OPTARG};;
		r) RUN=true;;
    m) MER=true;;
    o) PRO=true;;
    p) PLT=true;;
    g) GIF=true;;
    v) VID=true;;
	esac
done

# check to see if arguments were passed
if [ -z "$NAME" ]
then
	NAME=$DATETIME
	echo "-n, No name specified, using NAME=$NAME"
else
  echo "-n, Name specified, using NAME=$NAME"
fi
if [ -z "$CORES" ]
then
	CORES=1
	echo "-c, No number of cores specified, using CORES=$CORES"
else
  echo "-c, Number of cores specified, using CORES=$CORES"
fi
if [ "$RUN" = true ]
then
	echo "-r, Run the simulation"
fi
if [ "$MER" = true ]
then
	echo "-m, Merge the data"
fi
if [ "$PRO" = true ]
then
	echo "-o, Post-process data"
fi
if [ "$PLT" = true ]
then
	echo "-p, Plot the simulation"
fi
if [ "$GIF" = true ]
then
	echo "-g, Create gif of simulation"
fi
if [ "$VID" = true ]
then
	echo "-v, Create video (mp4) of the simulation"
fi
