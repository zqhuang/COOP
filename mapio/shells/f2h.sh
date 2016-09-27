mapdir=actlatest
outdir=actpost
maskfile=mask.fits
fwhm=${3}
lmin=${4}
lmax=2500
beam=${5}
lxcut=${6}
lycut=${7}
overwrite=${8}

#
./FSmooth -map ${mapdir}/${1}_I.fits -out ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}_I.fsm  -beam ${mapdir}/${beam} -fwhm ${fwhm} -lmin ${lmin} -lmax ${lmax} -mask ${mapdir}/${maskfile} -field I -lxcut ${lxcut} -lycut ${lycut} -overwrite ${overwrite}
./FSmoothQU -qmap ${mapdir}/${1}_Q.fits -umap ${mapdir}/${1}_U.fits -out_prefix ${outdir}/${2}_${fwhm}a_l${lmin}-${lmax}  -beam ${mapdir}/${beam} -fwhm ${fwhm} -lmin ${lmin} -lmax ${lmax}  -mask ${mapdir}/${maskfile} -field Q  -lxcut ${lxcut} -lycut ${lycut} -smoothmask act_matcoadd_standard_weight.fits  -overwrite ${overwrite}
