#!/bin/bash

# Run this script to submit a job to Graham
# Takes in optional arguments:
#	$ sh HPC_job_submit.sh -e <name of experiment> 				 Default: current datetime
#								  			 -s <number of simulations>      Default: 5
#								  			 -c <number of cores>            Default: 1
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

###############################################################################
# Miscellaneous variables

# Name of the job to submit
JOBNAME=$EXP

# Pull the most recent changes from git
git pull

# Prepare experiment directory
bash prepare_exp.sh -e $EXP -s $SIMS

###############################################################################
# Submit job to queue
# if [ "$ARGS" = true ]
# then
# 	sbatch --job-name=$JOBNAME $LANCEUR -e $EXP -s $SIMS -c $CORES -$RUN$MER$PRO$PLT$GIF$VID
# else
# 	sbatch --job-name=$JOBNAME $LANCEUR -e $EXP -s $SIMS -c $CORES
# fi
