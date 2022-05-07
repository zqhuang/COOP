#! /bin/bash
mapdir=actpost
peakprefix=${1}
prefix=${2}
fwhm=${3}
lmin=${4}
lmax=2500
radius=1.5
res=40 
nu=0
npw=0
norm=F
outprefix=n${norm}${npw}_${prefix}_on_${peakprefix}
outdir=actpics
output=${outdir}/${outprefix}_l${lmin}_${fwhm}a


###do theory
###get peaks
./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_I.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -peak RANDOM -orient RANDOM -nu ${nu} -hot T -cold T -maxn ${11} 

#./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nu0.dat -peak RANDOM -orient RANDOM -nu ${nu}  -hot T -cold T  -maxn ${11} 

#./FGetPeaks -map ${mapdir}/${peakprefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -out spots/${outprefix}_l${lmin}_${fwhm}a_Bhotcold_nu0.dat -peak RANDOM -orient RANDOM -nu ${nu} -hot T -cold T  -maxn ${11}


## stack with random rotation
### T on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_I.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -field T -out ${output}_randrot_TonThotcoldnu0 -width 4.5 -height 4 -min ${5} -max ${6}  -radius ${radius} -randrot T -res ${res}  -norm_to_corr ${norm} -norm_power ${npw}

#../utils/fasy.sh ${output}_randrot_TonThotcoldnu0.txt


###E on T
./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -field E -out ${output}_randrot_EonThotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}  -norm_to_corr ${norm} -norm_power ${npw}


####Q_r on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_QU.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -field QrUr -out ${output}_randrot_QrUronThotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}  -norm_to_corr ${norm} -norm_power ${npw}


#../utils/fasy.sh ${output}_randrot_EonThotcoldnu0.txt

### E on E
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_E.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nu0.dat -field E -out ${output}_randrot_EonEhotcoldnu0 -width 4.5 -height 4 -min ${9} -max ${10} -radius ${radius} -randrot T -res ${res}  -norm_to_corr ${norm} -norm_power ${npw}

#../utils/fasy.sh ${output}_randrot_EonEhotcoldnu0.txt


###B on T
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -field B -out ${output}_randrot_BonThotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res} -norm_to_corr ${norm} -norm_power ${npw} -want_caption F

#../utils/fasy.sh ${output}_randrot_BonThotcoldnu0.txt

###B on E
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Ehotcold_nu0.dat -field B -out ${output}_randrot_BonEhotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res}  -want_caption F -want_arrow F -norm_to_corr ${norm} -norm_power ${npw}
#../utils/fasy.sh ${output}_randrot_BonEhotcoldnu0.txt

##B on B
#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_B.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Bhotcold_nu0.dat -field B -out ${output}_randrot_BonBhotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -radius ${radius} -randrot T -res ${res} -norm_to_corr ${norm} -norm_power ${npw}
#../utils/fasy.sh ${output}_randrot_BonBhotcoldnu0.txt


##QrUr on T

#./FStack -map ${mapdir}/${prefix}_${fwhm}a_l${lmin}-${lmax}_QU.fsm -peaks spots/${outprefix}_l${lmin}_${fwhm}a_Thotcold_nu0.dat -field QrUr -out ${output}_randrot_QrUronThotcoldnu0 -width 4.5 -height 4 -min ${7} -max ${8} -min2 ${7} -max2 ${8} -radius ${radius} -randrot T -res ${res}  -norm_to_corr ${norm} -norm_power ${npw}

#mv ${output}_randrot_QrUronThotcoldnu0_1.txt ${output}_randrot_QronThotcoldnu0.txt
#mv ${output}_randrot_QrUronThotcoldnu0_2.txt ${output}_randrot_UronThotcoldnu0.txt

#../utils/fasy.sh ${output}_randrot_QronThotcoldnu0.txt
#../utils/fasy.sh ${output}_randrot_UronThotcoldnu0.txt



