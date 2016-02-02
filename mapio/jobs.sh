l1=250
l2=500
l3=1000
fwhm1=5
fwhm2=5
fwhm3=3
maxn1=25000
maxn2=100000
maxn3=400000
tmin=-15
tmax=30
qmin=-0.3
qmax=0.35
emin=-0.8
emax=3.2

./f2h.sh planck_smica planck ${fwhm1} ${l1} planck_beam.txt 0 0
./f2h.sh planck_smica planck ${fwhm2} ${l2} planck_beam.txt 0 0
./f2h.sh planck_smica planck ${fwhm3} ${l3} planck_beam.txt 0 0
./f2h.sh planck_smica xycut_planck ${fwhm1} ${l1} planck_beam.txt 50 80

./f2h.sh deep56_coadd act ${fwhm1} ${l1} beam_7ar2.txt 0 0
./f2h.sh deep56_coadd act ${fwhm2} ${l2} beam_7ar2.txt 0 0
./f2h.sh deep56_coadd act ${fwhm3} ${l3} beam_7ar2.txt 0 0
./f2h.sh deep56_coadd xycut_act ${fwhm1} ${l1} beam_7ar2.txt 50 80

./f2h.sh sim_1 noiseless ${fwhm1} ${l1} beam_7ar2.txt 0 0
./f2h.sh sim_1 noiseless  ${fwhm2} ${l2} beam_7ar2.txt 0 0
./f2h.sh sim_1 noiseless  ${fwhm3} ${l3} beam_7ar2.txt 0 0
./f2h.sh sim_1 xycut_noiseless ${fwhm1} ${l1} beam_7ar2.txt 50 80

./f2h.sh sim_with_noise_1 sim1 ${fwhm1} ${l1} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_1 sim1  ${fwhm2} ${l2} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_1 sim1  ${fwhm3} ${l3} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_1 xycut_sim1 ${fwhm1} ${l1} beam_7ar2.txt 50 80

./f2h.sh sim_with_noise_2 sim2 ${fwhm1} ${l1} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_2 sim2  ${fwhm2} ${l2} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_2 sim2  ${fwhm3} ${l3} beam_7ar2.txt 0 0
./f2h.sh sim_with_noise_2 xycut_sim2 ${fwhm1} ${l1} beam_7ar2.txt 50 80


./fullrun.sh planck act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh planck act ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./fullrun.sh xycut_planck xycut_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh planck act ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}

./fullrun.sh act act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh act act ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./fullrun.sh xycut_act xycut_act ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh act act ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}

./fullrun.sh sim1 sim1 ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh sim1 sim1 ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./fullrun.sh xycut_sim1 xycut_sim1 ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh sim1 sim1 ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}


./fullrun.sh noiseless noiseless ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh noiseless noiseless ${fwhm2} ${l2} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn2}
./fullrun.sh xycut_noiseless xycut_noiseless ${fwhm1} ${l1} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn1}
./fullrun.sh noiseless noiseless ${fwhm3} ${l3} ${tmin} ${tmax} ${qmin} ${qmax} ${emin} ${emax} ${maxn3}
