#!/bin/bash

# Run this script to submit a job to Graham
# Takes in optional arguments:
#	$ sh HPC_job_submit.sh -e <name of experiment> 				 Default: current datetime
#								  			 -s <number of simulations>      Default: 5
#								  			 -c <number of cores>            Default: 32
#								  			 -r <run simulation>             Default: False
#             	  			 -m <merge h5 files>             Default: False
#								  			 -o <post-process data>          Default: False
#								  			 -p <plot simulation>            Default: False
#             	  			 -g <create gif>                 Default: False
#             	  			 -v <create video, mp4>          Default: False

# Current datetime
DATETIME=`date +"%Y-%m-%d_%Hh%M"`

# Having a ":" after a flag means an option is required to invoke that flag
while getopts e:s:c:rmopgv option
do
	case "${option}"
		in
		e) EXP=${OPTARG};;
		s) SIMS=${OPTARG};;
		c) CORES=${OPTARG};;
		r) RUN=r;;
    m) MER=m;;
    o) PRO=o;;
    p) PLT=p;;
    g) GIF=g;;
    v) VID=v;;
	esac
done

# check to see if arguments were passed
if [ -z "$EXP" ]
then
	EXP=$DATETIME
	echo "-e, No name specified, using EXP=$EXP"
else
  echo "-e, Name specified, using EXP=$EXP"
fi
if [ -z "$SIMS" ]
then
	SIMS=5
	echo "-s, No number of simulations specified, using SIMS=$SIMS"
else
  echo "-s, Number of simulations specified, using SIMS=$SIMS"
fi
if [ -z "$CORES" ]
then
	CORES=32
	echo "-c, No number of cores specified, using CORES=$CORES"
else
  echo "-c, Number of cores specified, using CORES=$CORES"
fi
if [ "$RUN" = r ]
then
	ARGS=true
	echo "-r, Run the simulation"
fi
if [ "$MER" = m ]
then
	ARGS=true
	echo "-m, Merge the data"
fi
if [ "$PRO" = o ]
then
	ARGS=true
	echo "-o, Post-process data"
fi
if [ "$PLT" = p ]
then
	ARGS=true
	echo "-p, Plot the simulation"
fi
if [ "$GIF" = g ]
then
	ARGS=true
	echo "-g, Create gif of simulation"
fi
if [ "$VID" = v ]
then
	ARGS=true
	echo "-v, Create video (mp4) of the simulation"
fi

###############################################################################
JOBNAME=$EXP
LANCEUR="HPC_lanceur.slrm"
###############################################################################
# Pull the most recent changes from git
git pull

# Check if directory for experiment exists
if [ ! -d _experiments/${JOBNAME} ]
then
	mkdir _experiments/${JOBNAME}
fi

# Submit job to queue
if [ "$ARGS" = true ]
then
	sbatch --job-name=$JOBNAME $LANCEUR -e $EXP -s $SIMS -c $CORES -$RUN$MER$PRO$PLT$GIF$VID
else
	sbatch --job-name=$JOBNAME $LANCEUR -e $EXP -s $SIMS -c $CORES
fi
