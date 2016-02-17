l1=250
l2=500
l3=1000
fwhm1=5
fwhm2=5
fwhm3=3
maxn1=90000
maxn2=100000
maxn3=400000
tmin=-15
tmax=30
qmin=-0.3
qmax=0.35
emin=-0.8
emax=3.2
lxcut=${1}
lycut=${2}

#./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm1} ${l1} planck_beam.txt ${lxcut} ${lycut}
./f2h.sh deep56_coadd cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}


#./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_planck ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
