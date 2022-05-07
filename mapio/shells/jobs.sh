fwhm=5
tmin=-5
tmax=15
qmin=-0.1
qmax=0.08
emin=-1.2
emax=3.2
lxcut=90
lycut=50

ell=350

maxn=250000
./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm} ${ell} planck_beam.txt ${lxcut} ${lycut} F
./f2h.sh act_S2_map cutx${lxcut}y${lycut}_act ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F

./FSim -isim 1
./submap_spn.sh 1
./f2h.sh sim_1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
./f2h.sh spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F

./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_s1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}

maxn=100000
for isim in `seq 2 30`
do
    ./FSim -isim ${isim}
    ./submap_spn.sh ${isim}
    ./f2h.sh sim_${isim} cutx${lxcut}y${lycut}_s${isim} ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
    ./f2h.sh spn${isim} cutx${lxcut}y${lycut}_spn${isim} ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
    ./partrun.sh cutx${lxcut}y${lycut}_s${isim} cutx${lxcut}y${lycut}_s${isim} ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
    ./partrun.sh cutx${lxcut}y${lycut}_spn${isim} cutx${lxcut}y${lycut}_spn${isim} ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
    rm -f actpost/*_s${isim}_* actpost/*_spn${isim}_*
done


ell=800

maxn=250000
./f2h.sh planck_smica cutx${lxcut}y${lycut}_planck ${fwhm} ${ell} planck_beam.txt ${lxcut} ${lycut} F
./f2h.sh act_S2_map cutx${lxcut}y${lycut}_act ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F

./FSim -isim 1
./submap_spn.sh 1
./f2h.sh sim_1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
./f2h.sh spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F

./partrun.sh cutx${lxcut}y${lycut}_planck cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_s1 cutx${lxcut}y${lycut}_s1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_act cutx${lxcut}y${lycut}_act ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
./partrun.sh cutx${lxcut}y${lycut}_spn1 cutx${lxcut}y${lycut}_spn1 ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}

maxn=100000
for isim in `seq 2 30`
do
    ./FSim -isim ${isim}
    ./submap_spn.sh ${isim}
    ./f2h.sh sim_${isim} cutx${lxcut}y${lycut}_s${isim} ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
    ./f2h.sh spn${isim} cutx${lxcut}y${lycut}_spn${isim} ${fwhm} ${ell} beam_gaussian_fwhm0arcmin.txt ${lxcut} ${lycut} F
    ./partrun.sh cutx${lxcut}y${lycut}_s${isim} cutx${lxcut}y${lycut}_s${isim} ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
    ./partrun.sh cutx${lxcut}y${lycut}_spn${isim} cutx${lxcut}y${lycut}_spn${isim} ${fwhm} ${ell} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn}
    rm -f actpost/*_s${isim}_* actpost/*_spn${isim}_*
done

