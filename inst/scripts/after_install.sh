#!/bin/bash
. $(dirname $0)/constants.sh
# move into the test folder
cd $source_path/$code_deploy_path/$test_folder
# # execute infra and config files
sudo chmod +x utils/install.R
# # to input '1' infinitely when the console prompts for input a
yes 1 | ./utils/install.R
# make the file executable
sudo chmod +x dataprocess.R
# & to run in background
# pointing output to dev/null so that code deploy doesn't stall
cmd_arg=$source_path/$code_deploy_path/$test_folder/$test_utils_folder
yes a | ./dataprocess.R ${cmd_arg} > /dev/null 2> /dev/null < /dev/null &