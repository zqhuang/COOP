#! /bin/bash
echo ${1} > asyplot.config
if [[ -n ${2} ]]
then
    output=${2}.pdf
else
  if [[ $1 =~ \.txt$ ]]
  then
      output=${1/\.txt/\.pdf}
  else
      echo "$1 does not have a postfix .txt"
      exit
  fi
fi
asy ${HOME}/work/GitHub/COOP/utils/asyplot.asy -o $output
rm -f asyplot.config

