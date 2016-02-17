#! /bin/bash
mapdir=actpost
peakprefix=${1}
prefix=${2}
fwhm=${3}
lmin=${4}
lmax=4000
radius=`echo 500./${lmin} | bc -l`
res=$(( 100000 / ${fwhm} / (${lmin}+ 150) ))
nu=0.5
outprefix=${prefix}_on_${peakprefix}
outdir=acttalk
output=${outdir}/${outprefix}_l${lmin}_${fwhm}a

###do theory
#./GetTheo -cl lcdm_lensedCls.dat -auto 1 -cross 1 -field T -peak RANDOM -out ${outdir}/theory_noiseless_l${lmin}_${fwhm}a_TonHotnupt5 -nu ${nu} -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-20)) -hpauto_highl $((${lmin}+20)) -hpcross_lowl $((${lmin}-20)) -hpcross_highl $((${lmin}+20)) -radius ${radius}

#./GetTheo -cl lcdm_lensedCls.dat -auto 1 -cross 4 -field E -peak RANDOM -out ${outdir}/theory_noiseless_l${lmin}_${fwhm}a_EonHotnupt5 -nu ${nu} -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-20)) -hpauto_highl $((${lmin}+20)) -hpcross_lowl $((${lmin}-20)) -hpcross_highl $((${lmin}+20)) -radius ${radius}

#./GetTheo -cl lcdm_lensedCls.dat -auto 2 -cross 2 -field E -peak RANDOM -out ${outdir}/theory_noiseless_l${lmin}_${fwhm}a_EonEhotnupt5 -nu ${nu} -orient RANDOM -fwhm ${fwhm} -hot T -lmax ${lmax} -hpauto_lowl $((${lmin}-20)) -hpauto_highl $((${lmin}+20)) -hpcross_lowl $((${lmin}-20)) -hpcross_highl $((${lmin}+20)) -radius ${radius}

###get peaks
./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_I.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nupt5.dat -peak RANDOM -orient RANDOM -nu ${nu} -hot T -cold T -maxn ${11}

#./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nupt5.dat -peak RANDOM -orient RANDOM -nu ${nu}  -hot T -cold T  -maxn ${11}

#./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Bhotcold_nupt5.dat -peak RANDOM -orient RANDOM -nu ${nu} -hot T -cold T  -maxn ${11}


## stack with random rotation
### T on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_I.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nupt5.dat -field T -out ${output}_randrot_TonThotcoldnupt5 -width 4.5 -height 4 -min ${5} -max ${6}  -radius ${radius} -randrot T -res ${res} 

#../utils/fasy.sh ${output}_randrot_TonThotcoldnupt5.txt


###E on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nupt5.dat -field E -out ${output}_randrot_EonThotcoldnupt5 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}

#../utils/fasy.sh ${output}_randrot_EonThotcoldnupt5.txt

### E on E
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nupt5.dat -field E -out ${output}_randrot_EonEhotcoldnupt5 -width 4.5 -height 4 -min ${9} -max ${10} -radius ${radius} -randrot T -res ${res}

#../utils/fasy.sh ${output}_randrot_EonEhotcoldnupt5.txt


###B on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nupt5.dat -field B -out ${output}_randrot_BonThotcoldnupt5 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}

#../utils/fasy.sh ${output}_randrot_BonThotcoldnupt5.txt

###B on E
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nupt5.dat -field B -out ${output}_randrot_BonEhotcoldnupt5 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}
#../utils/fasy.sh ${output}_randrot_BonEhotcoldnupt5.txt

##B on B
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Bhotcold_nupt5.dat -field B -out ${output}_randrot_BonBhotcoldnupt5 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}
#../utils/fasy.sh ${output}_randrot_BonBhotcoldnupt5.txt


##QrUr on T

./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_QU.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nupt5.dat -field QrUr -out ${output}_randrot_QrUronThotcoldnupt5 -width 4.5 -height 4 -min ${7} -max ${8} -min2 ${7} -max2 ${8} -radius ${radius} -randrot T -res ${res}

mv ${output}_randrot_QrUronThotcoldnupt5_1.txt ${output}_randrot_QronThotcoldnupt5.txt
mv ${output}_randrot_QrUronThotcoldnupt5_2.txt ${output}_randrot_UronThotcoldnupt5.txt

../utils/fasy.sh ${output}_randrot_QronThotcoldnupt5.txt
../utils/fasy.sh ${output}_randrot_UronThotcoldnupt5.txt



