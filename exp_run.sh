#!/bin/bash
# A bash script to run many simulations in an experiment
# Takes in optional arguments:
#	$ sh exp_run.sh -e <name of experiment> 				Default: current datetime
#								  -s <number of simulations>      Default: 1
#								  -c <number of cores>            Default: 1
#									-a <ask before deleting>				Default: False
#								  -r <run simulation>             Default: False
#             	  -m <merge h5 files>             Default: False
#								  -o <post-process data>          Default: False
#								  -p <plot simulation>            Default: False
#             	  -g <create gif>                 Default: False
#             	  -v <create video, mp4>          Default: False

# Current datetime
DATETIME=`date +"%Y-%m-%d_%Hh%M"`

# Having a ":" after a flag means an option is required to invoke that flag
while getopts e:s:c:armopgv option
do
	case "${option}"
		in
		e) EXP=${OPTARG};;
		s) SIMS=${OPTARG};;
		c) CORES=${OPTARG};;
		a) ASK=true;;
		r) RUN=r;;
    m) MER=m;;
    o) PRO=o;;
    p) PLT=p;;
    g) GIF=g;;
    v) VID=v;;
	esac
done

ARGS=false
# check to see what arguments were passed
if [ -z "$EXP" ]
then
	EXP="test_exp"
	echo "-e, No name specified, using EXP=$EXP"
else
  echo "-e, Name specified, using EXP=$EXP"
fi
if [ -z "$SIMS" ]
then
	SIMS=1
	echo "-s, No number of simulations specified, using SIMS=$SIMS"
else
  echo "-s, Number of simulations specified, using SIMS=$SIMS"
fi
if [ -z "$CORES" ]
then
	CORES=1
	echo "-c, No number of cores specified, using CORES=$CORES"
else
  echo "-c, Number of cores specified, using CORES=$CORES"
fi
if [ "$ASK" = a ]
then
	ARGS=true
	echo "-a, Will ask before overwriting files"
else
	echo "-a, Will not ask before overwriting files"
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
# Path to experiment directory
EXP_PATH="_experiments/${EXP}"

# Prepare experiment directory
bash prepare_exp.sh -e $EXP -s $SIMS

###############################################################################
# Run simulations

# Navigate to experiment directory
cd ${EXP_PATH}
# Loop over all simulations if any of -rm args are activated
if [ "$RUN" = r ] || [ "$MER" = m ]
then
	for (( id=0; id<$SIMS; id++ ))
	do
		echo "--Running simulation $id--"
		# Format ID with correct number of digits
		printf -v PID "%03d" $id
		# Navigate to the simulation directory
		cd _simulations/${PID}_$EXP
		# Check whether any switch arguments were activated
		if [ "$ARGS" = true ]
		then
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES -$ASK$RUN$MER$PRO$PLT$GIF$VID
		else
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES
		fi
		# Back out to experiment directory
		cd ..
		cd ..
	done
fi

# exit 0
###############################################################################
# The version of python to use
python_command="python3"
# Name of csv data file
csv_data_file="sim_data.csv"
# Name of plotting script
plot_data_file="plot_exp_data.py"
###############################################################################
# Plot experiments' transmission coefficients if simulations were run
if [ "$PRO" = o ]
then
	echo ''
	echo '--Plotting transmission coefficients--'
	# Check to make sure snapshots exists
	echo "Checking for csv data file"
	if [ -e _simulations/000_${EXP}/$csv_data_file ]
	then
		echo "Data found"
	else
		echo "Cannot find data. Aborting script"
		exit 1
	fi
	echo 'Plotting experiment data'
	${python_command} $plot_data_file $EXP $SIMS
	echo 'Done plotting'
fi


echo ''
echo 'Done running experiment'
echo ''
