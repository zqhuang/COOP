#! /bin/bash
mapdir=act16
prefix=${1}
lmin=${2}
lmax=2500
fwhm=5
radius=1.5
outprefix=${prefix}

output=ACTstacking/${outprefix}_${fwhm}a_l${lmin}


#./GetPeaks -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -out spots/${outprefix}_l${lmin}_${fwhm}a_hot_nu1.dat -peak RANDOM -orient RANDOM -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -nu 1



./Stack -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${outprefix}_l${lmin}_${fwhm}a_hot_nu1.dat -field T -out ${output}_randrot_TonHotnu1 -width 4.5 -height 4 -min ${3} -max ${4} -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius ${radius} -randrot T -res 50

../utils/fasy.sh ${output}_randrot_TonHotnu1.txt


./Stack -map ${mapdir}/${prefix}_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${outprefix}_l${lmin}_${fwhm}a_hot_nu1.dat -field T -out ${output}_norot_TonHotnu1 -width 4.5 -height 4 -min ${3} -max ${4} -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius ${radius} -randrot F -res 50

../utils/fasy.sh ${output}_norot_TonHotnu1.txt


./Stack -map ${mapdir}/${prefix}_qu_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${outprefix}_l${lmin}_${fwhm}a_hot_nu1.dat -field QrUr -out ${output}_randrot_QrUronHotnu1 -width 4.5 -height 4 -min ${5} -max ${6}   -min2 ${5} -max2 ${6}  -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius ${radius} -randrot T -res 50

../utils/fasy.sh ${output}_randrot_QrUronHotnu1_1.txt
../utils/fasy.sh ${output}_randrot_QrUronHotnu1_2.txt
mv ${output}_randrot_QrUronHotnu1_1.pdf  ${output}_randrot_QronHotnu1.pdf
mv ${output}_randrot_QrUronHotnu1_2.pdf  ${output}_randrot_UronHotnu1.pdf

./Stack -map ${mapdir}/${prefix}_qu_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${outprefix}_l${lmin}_${fwhm}a_hot_nu1.dat -field QrUr -out ${output}_norot_QrUronHotnu1 -width 4.5 -height 4  -min ${5} -max ${6}  -min2 ${5} -max2 ${6}  -mask ${mapdir}/${prefix}_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius ${radius} -randrot F -res 50

../utils/fasy.sh ${output}_norot_QrUronHotnu1_1.txt
../utils/fasy.sh ${output}_norot_QrUronHotnu1_2.txt
mv ${output}_norot_QrUronHotnu1_1.pdf  ${output}_norot_QronHotnu1.pdf
mv ${output}_norot_QrUronHotnu1_2.pdf  ${output}_norot_UronHotnu1.pdf


#./GetTheo -cl planck14best_lensedCls.dat -auto 1 -cross 1 -field T -peak RANDOM -out actstack_zqcoadd/theory_noiseless_l${lmin}_${fwhm}a_TonHotnu1 -nu 1 -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-10)) -hpauto_highl $((${lmin}+10)) -hpcross_lowl $((${lmin}-10)) -hpcross_highl $((${lmin}+10)) -radius ${radius}

#./GetTheo -cl planck14best_lensedCls.dat -auto 1 -cross 4 -field E -peak RANDOM -out actstack_zqcoadd/theory_noiseless_l${lmin}_${fwhm}a_EonHotnu1 -nu 1 -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-10)) -hpauto_highl $((${lmin}+10)) -hpcross_lowl $((${lmin}-10)) -hpcross_highl $((${lmin}+10)) -radius ${radius}

