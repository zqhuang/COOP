#!/bin/bash
if [ ${1} == 'f90' ] 
then
    for i in `ls *.f90`
      do 
      echo ${i}
      cat -n ${i} |grep -i --ignore-case ${2}
    done
    for i in `ls *.F90`
      do 
      echo ${i}
      cat -n ${i} |grep -i --ignore-case ${2}
    done
else
    for i in `ls *.${1}`
      do 
      echo ${i}
      cat -n ${i} |grep -i --ignore-case ${2}
    done
fi

