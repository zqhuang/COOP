l1=250
l2=500
l3=1000
fwhm1=5
fwhm2=5
fwhm3=3
maxn1=50000
maxn2=100000
maxn3=250000
tmin=-15
tmax=30
qmin=-0.3
qmax=0.3
emin=-0.8
emax=3.2
lxcut=90
lycut=50


#./f2h.sh deep56_coadd cutx${lxcut}y${lycut}_act ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}
./f2h.sh noise1 cutx${lxcut}y${lycut}_nb1 ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}


#./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
#./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
./partrun.sh cutx${lxcut}y${lycut}_nb1 cutx${lxcut}y${lycut}_nb1 ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
