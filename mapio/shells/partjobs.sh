ell=800
fwhm=5
maxn=800000
tmin=-15
tmax=30
qmin=-0.1
qmax=0.08
emin=-1.2
emax=3.2
lxcut=90
lycut=50


./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm} ${ell} planck_beam.txt ${lxcut} ${lycut} F
./f2h.sh act_S2_map cutx${lxcut}y${lycut}_act ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
./f2h.sh sim_1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
./f2h.sh spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F

./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_s1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}

