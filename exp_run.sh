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
		a) ASK=a;;
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
# The command and arguments for running scripts with mpi
mpiexec_command="mpiexec"
# The version of python to use
python_command="python3"
# Path to snapshot files
snapshot_name="snapshots"
snapshot_path="${NAME}/${snapshot_name}"
# Name of output directory
output_dir="outputs"
# Path to frames
frames_path='frames'
# Clean up the snapshots after merging
cleanup_snapshots="True"

# Name of the script to write the parameters
params_script='write_params.py'
params_file='params.py'
# Name of the main code file
code_file='main.py'
# Name of switchboard file
switchboard="switchboard"
switch_file="${switchboard}.py"
sim_switch="${NAME}/${NAME}_${switch_file}"
# Name of merging file
merge_file="merge.py"
# Helper code files
helper_funcs="helper_functions.py"
helper_funcs_CD="helper_functions_CD.py"
# Name of post processing file
post_process="post_process.py"
# Name of experiment plotting file
plot_exp="plot_exp_data.py"
# Name of slice plotting file
plot_file="plot_slices.py"
# Name of gif creation file
gif_cre_file="create_gif.py"
# Group all the code files for ease of calling
CODE_FILES="$params_script $code_file $switch_file $merge_file $helper_funcs $helper_funcs_CD $post_process $plot_exp $plot_file $gif_cre_file"

###############################################################################
echo ''
echo '--Checking experiment directory--'
# Check if experiments folder exists
if [ -e _experiments ]
then
	echo 'Experiment folder exists'
else
	echo 'Experiment folder not found. Aborting script'
	exit 1
fi
echo ''
OVERWRITE_CODE_FILES=true
# Check if this experiment has been created
if [ -e _experiments/$EXP ]
then
	echo "Experiment for $EXP exists"
	if [ "$ASK" = false ]
	then
		echo "Overwriting code files"
		# Go in to directory, remove code files, come back out to main directory
		cd _experiments/$EXP
		rm -rf $CODE_FILES
		cd ..
		cd ..
	else
		# Begin loop waiting for user to confirm
		while true
		do
			read -r -p "Overwrite old code files in $EXP [y/n]? Or cancel? [Ctrl+c] " input
			case $input in
				[yY][eE][sS]|[yY])
			echo "Yes"
			# Go in to directory, remove code files, come back out to main directory
			cd _experiments/$EXP
			rm -rf $CODE_FILES
			cd ..
			cd ..
			break
			;;
				[nN][oO]|[nN])
			echo "No"
			OVERWRITE_CODE_FILES=false
			break
					;;
				*)
			echo "Invalid input"
			;;
			esac
		done
	fi
else
	echo "Creating experiment for $EXP"
	mkdir "_experiments/$EXP"
fi
if [ $OVERWRITE_CODE_FILES = true ] || [ "$ASK" = false ]
then
	echo 'Copying code files to experiment directory'
	cp $CODE_FILES _experiments/${EXP}
fi
echo ''

###############################################################################
###############################################################################
# run simulations
for (( id=0; id<$SIMS; id++ ))
do
	echo "--Running simulation $id--"
	# If -c is more than 1, it will run the simulations in parallel
	if [ $CORES -eq 1 ] # run serially
	then
		# Check whether any switch arguments were activated
		if [ "$ARGS" = true ]
		then
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES -$ASK$RUN$MER$PRO$PLT$GIF$VID
		else
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES
		fi
	else # run in parallel
		# Check whether any switch arguments were activated
		#		The `&` will mean all simulations will start at around the same time
		if [ "$ARGS" = true ]
		then
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES -$ASK$RUN$MER$PRO$PLT$GIF$VID &
		else
			bash run.sh -e $EXP -i $id -s $SIMS -c $CORES &
		fi
	# Record pid of last command run
	pids[${id}]=$!
	fi
done

# Wait for all pids to complete before continuing
for pid in ${pids[*]}; do
    wait $pid
done

# exit 0
###############################################################################
# The version of python to use
python_command="python3"
# Name of csv data file
csv_data_file="exp_data.csv"
# Name of plotting script
plot_data_file="plot_exp_data.py"
# Name of switchboard file
switchboard="switchboard"
###############################################################################
# Plot experiments' transmission coefficients if simulations were run
if [ "$PRO" = o ]
then
	cd _experiments/${EXP}
	echo ''
	echo '--Plotting transmission coefficients--'
	# Check to make sure snapshots exists
	echo "Checking for csv data file"
	if [ -e $csv_data_file ]
	then
		echo "Data found"
	else
		echo "Cannot find data. Aborting script"
		exit 1
	fi
	echo 'Plotting experiment data'
	${python_command} $plot_data_file $EXP $switchboard
	echo 'Done plotting'
fi


echo ''
echo 'Done running experiment'
echo ''
