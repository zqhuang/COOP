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
qmin=-0.4
qmax=0.4
emin=-0.8
emax=3.2
lxcut=90
lycut=50

#./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh act_matcoadd cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}

#./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
#./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
#./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}


#./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm2} ${l2} planck_beam.txt ${lxcut} ${lycut}
#./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh act_matcoadd cutx${lxcut}y${lycut}_act ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}

#./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
#./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
#./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}

./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm2} ${l2} planck_beam.txt ${lxcut} ${lycut}
./f2h.sh act_matcoadd cutx${lxcut}y${lycut}_act ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}
./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}
./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm2} ${l2} beam_7ar2.txt ${lxcut} ${lycut}

./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_act ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
