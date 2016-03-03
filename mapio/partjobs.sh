l1=250
l2=500
l3=1000
fwhm1=5
fwhm2=5
fwhm3=3
maxn1=80000
maxn2=100000
maxn3=250000
tmin=-15
tmax=30
qmin=-0.4
qmax=0.4
emin=-0.8
emax=3.2
lxcut=0
lycut=0

#./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh act_matcoadd cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}

#./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
#./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
#./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}


#./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm3} ${l3} planck_beam.txt ${lxcut} ${lycut}
#./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh act_matcoadd cutx${lxcut}y${lycut}_act ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}
#./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm3} ${l3} beam_7ar2.txt ${lxcut} ${lycut}

#./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
#./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
#./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}

#./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm1} ${l1} planck_beam.txt ${lxcut} ${lycut}
#./f2h.sh deep56nopf cutx${lxcut}y${lycut}_actnopf ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut} T
#./f2h.sh noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut} F
./f2h.sh sim_1 cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} beam_7ar2.txt ${lxcut} ${lycut}  T

#./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./partrun.sh cutx${lxcut}y${lycut}_noiseless cutx${lxcut}y${lycut}_noiseless ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
#./partrun.sh cutx${lxcut}y${lycut}_actnopf cutx${lxcut}y${lycut}_actnopf ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1} 
#./partrun.sh cutx${lxcut}y${lycut}_noise1 cutx${lxcut}y${lycut}_noise1 ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
