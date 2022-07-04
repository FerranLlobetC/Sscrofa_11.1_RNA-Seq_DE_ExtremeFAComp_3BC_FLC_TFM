#!/bin/bash

: " Input parameters "
# File that requires it's shebang to be updated:
FILE=$1
# Programing language name (bin) to include a shebang of:
LANG=$2

: " Safe-lock to ensure that only '.py' files are modified "

if [[ $FILE != *.py ]]; then
  echo -e "\033[1;48;2;103;7;78mOnly .py files can be processed by this script!!\033[0m"
  exit 1
fi

: " Saving of the program location directory "
# Addition of an error message when the file is not found:
if [[ $(command -v $LANG) ]]; then
  LNG_BIN=$(command -v $LANG)
else
  echo -e "\033[1;48;2;103;7;78mThe program is not in this machine!!\033[0m"
  exit 1
fi

: " Replacing the file shebang with the newly stored Python3 directio"
# Use of the powerful 'gsub' function of 'awk' to replace the shebang
awk -v Language=$LNG_BIN '{gsub("#![.-z]+|#!+", "#!"Language); print $0}' $FILE > $FILE.tmp

: " Update of the modified file "
rm $FILE
mv $FILE.tmp $FILE

: " Execution permission to the modified file "
chmod u+x $FILE

echo -e "\033[1;38;2;177;225;175m$LANG shebang of file $FILE successfully modified to match the current machine!!\033[0m \U0001F978"