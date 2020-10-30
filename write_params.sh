#!/bin/bash
# A bash script to write parameters to a python file
# Takes in optional arguments:
#	$ sh run.sh -x <param1, mL>         				Default: 0

# Having a ":" after a flag means an option is required to invoke that flag
while getopts x: option
do
	case "${option}"
		in
		x) P1=${OPTARG};;
	esac
done

# check to see what arguments were passed
if [ -z "$P1" ]
then
	P1=0
	echo "Parameter 1 not set, using $P1"
fi

###############################################################################
# The file to write to
params_file="parameters.py"

# Writing out the content between lines with EOL to the file
cat >> $params_file <<EOL
mL = $P1
theta = 0.7
EOL
