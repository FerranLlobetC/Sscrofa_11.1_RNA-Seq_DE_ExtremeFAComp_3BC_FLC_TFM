#!/bin/bash

: " Input parameters "
# File that requires it's shebang to be updated:
FILE=$1

: " Safe-lock to ensure that only '.py' files are modified "

if [[ $FILE != *.py ]]; then
  echo -e "\033[1;38;2;103;7;78mOnly .py files can be processed by this script!!\033[0m"
  exit 1
fi

: " Saving of the program location directory "
PY3=$(which python3)

: " Replacing the file shebang with the newly stored Python3 directio"
# Use of the powerful 'gsub' function of 'awk' to replace the shebang
awk -v Python3=$PY3 '{gsub("#![.-z]+", "#!"Python3); print $0}' $FILE > $FILE.tmp

: " Update of the modified file "
rm $FILE
mv $FILE.tmp $FILE

echo -e "Python3 shebang of file $FILE successfully modified to match the current machine!! \U0001F978"