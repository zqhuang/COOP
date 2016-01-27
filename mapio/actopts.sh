#! /bin/bash
mapdir=act16
lmin=${1}
lmax=2500
fwhm=3

output=actstack_zqcoadd/${mapdir}_l${lmin}_randrot_TonHotnu1

./GetPeaks -map ${mapdir}/act_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -out spots/${mapdir}l${lmin}_hot_nu1.dat -peak RANDOM -orient RANDOM -mask ${mapdir}/act_imask_${fwhm}a_l${lmin}-${lmax}.fits -nu 1

./Stack -map ${mapdir}/act_TQTUT_${fwhm}a_l${lmin}-${lmax}.fits -peaks spots/${mapdir}l${lmin}_hot_nu1.dat -field T -out ${output} -width 4.5 -height 4 -min ${2} -max ${3} -mask ${mapdir}/act_imask_${fwhm}a_l${lmin}-${lmax}.fits -radius 1 -randrot T -res 50

../utils/fasy.sh ${output}.txt
./GetTheo -cl planck14best_lensedCls.dat -auto 1 -cross 1 -field T -peak RANDOM -out actstack_zqcoadd/theory_noiseless_l${lmin}_TonHotnu1 -nu 1 -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-10)) -hpauto_highl $((${lmin}+10)) -hpcross_lowl $((${lmin}-10)) -hpcross_highl $((${lmin}+10)) -radius 1