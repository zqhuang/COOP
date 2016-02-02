mapdir=act16
outdir=actpost
fwhm=${3}
lmin=${4}
lmax=4000
beam=${5}
lxcut=${6}
lycut=${7}
#planck_beam.txt
#
./FSmooth -map ${mapdir}/${1}_I.fits -out ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_I.fsm  -beam ${mapdir}/${beam} -fwhm ${fwhm} -lmin ${lmin} -lmax ${lmax} -mask ${mapdir}/deep56_weight_I.fits -field I -lxcut ${lxcut} -lycut ${lycut}
./FSmooth -map ${mapdir}/${1}_Q.fits -out ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_Q.fsm  -beam ${mapdir}/${beam} -fwhm ${fwhm} -lmin ${lmin} -lmax ${lmax}  -mask ${mapdir}/deep56_weight_Q.fits -field Q  -lxcut ${lxcut} -lycut ${lycut}
./FSmooth -map ${mapdir}/${1}_U.fits -out ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_U.fsm  -beam ${mapdir}/${beam} -fwhm ${fwhm} -lmin ${lmin} -lmax ${lmax}  -mask ${mapdir}/deep56_weight_U.fits -field U  -lxcut ${lxcut} -lycut ${lycut}
./FMerge  -map1 ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_Q.fsm -map2 ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_U.fsm -out ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_QU.fsm
./FQU2EB  ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_QU.fsm  ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}
