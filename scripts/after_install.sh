#!/bin/bash
. $(dirname $0)/constants.sh
#Move into the test folder
cd $source_path/$code_deploy_path/$test_folder
# execute infra and config files
sudo chmod +x utils/infra.R
# to input '1' infinitely when the console prompts for input a
yes 1 | ./utils/infra.R
# make the file executable
sudo chmod +x dataprocess.R
#& to run in background
yes a | ./dataprocess.R > /dev/null 2> /dev/null < /dev/null &