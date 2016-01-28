#! /bin/bash
mapdir=act16
lmin=${1}
lmax=2500
fwhm=3
prefix=act

output=actstack_zqcoadd/${mapdir}_${fwhm}a_l${lmin}

./GetPeaks -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -out spots/${mapdir}l${lmin}_hot_nu1.dat -peak RANDOM -orient RANDOM -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -nu 1

./Stack -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${mapdir}l${lmin}_hot_nu1.dat -field T -out ${output}_randrot_TonHotnu1 -width 4.5 -height 4 -min ${2} -max ${3} -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius 1 -randrot T -res 50

../utils/fasy.sh ${output}_randrot_TonHotnu1.txt

./Stack -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${mapdir}l${lmin}_hot_nu1.dat -field T -out ${output}_norot_TonHotnu1 -width 4.5 -height 4 -min ${2} -max ${3} -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius 1 -randrot F -res 50

../utils/fasy.sh ${output}_norot_TonHotnu1.txt


./Stack -map ${mapdir}/${prefix}_qu_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${mapdir}l${lmin}_hot_nu1.dat -field QrUr -out ${output}_randrot_QrUronHotnu1 -width 4.5 -height 4  -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius 1 -randrot T -res 50

../utils/fasy.sh ${output}_randrot_QrUronHotnu1.txt


./Stack -map ${mapdir}/${prefix}_qu_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${mapdir}l${lmin}_hot_nu1.dat -field QrUr -out ${output}_norot_QrUronHotnu1 -width 4.5 -height 4  -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius 1 -randrot F -res 50

../utils/fasy.sh ${output}_norot_QrUronHotnu1.txt

#./GetTheo -cl planck14best_lensedCls.dat -auto 1 -cross 1 -field T -peak RANDOM -out actstack_zqcoadd/theory_noiseless_l${lmin}_TonHotnu1 -nu 1 -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-10)) -hpauto_highl $((${lmin}+10)) -hpcross_lowl $((${lmin}-10)) -hpcross_highl $((${lmin}+10)) -radius 1