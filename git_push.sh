#!/usr/bin/env bash
# Author: Mikhail Schee

# A script to automate git add, commit, and push

DATE=`date +"%y/%m/%d %T"`
COMMIT_M="\"$DATE update\""

set -x # echos each command as it is executed

git add --all
# eval required to avoid trouble with quotation marks
eval "git commit -m $COMMIT_M"
git push
