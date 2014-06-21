#! /bin/bash
rm -f tmp.eps
echo ${1} > asyplot.config
if [[ -n ${2} ]]
then
    epsfile=${2}.eps
else
  if [[ $1 =~ \.txt$ ]]
  then
    epsfile=${1/\.txt/\.eps}
  else
    epsfile=${1}\.eps
  fi
fi
asy ${HOME}/work/CosmoLib/utils/asyplot.asy -o tmp.eps
eps2eps tmp.eps $epsfile
epstopdf $epsfile
rm -f $epsfile tmp.eps
rm -f asyplot.config
