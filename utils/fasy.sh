#! /bin/bash
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
rm -f tmp_${epsfile}
asy ${HOME}/work/GitHub/COOP/utils/asyplot.asy -o tmp_${epsfile}
eps2eps tmp_${epsfile} ${epsfile}
epstopdf $epsfile
rm -f tmp_${epsfile} ${epsfile}
rm -f asyplot.config
