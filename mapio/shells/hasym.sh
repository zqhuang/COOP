#!/bin/bash
if [ -n "" ]
then   
./GetPeaks -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits  -mask planck14/HemAsym_south_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot T -out peaks/sasym_hot_nu0p5.dat
./GetPeaks -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -mask planck14/HemAsym_north_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot T -out peaks/nasym_hot_nu0p5.dat
./Stack -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -peaks peaks/sasym_hot_nu0p5.dat -field T -mask planck14/HemAsym_south_int_mask_120a_0128.fits -randrot T -out stacked/sasym_120a_hot_nu0p5
./Stack -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -peaks peaks/nasym_hot_nu0p5.dat -field T -mask planck14/HemAsym_north_int_mask_120a_0128.fits -randrot T -out stacked/nasym_120a_hot_nu0p5
./GetPeaks -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits  -mask planck14/HemAsym_south_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot F -out peaks/sasym_cold_nu0p5.dat
./GetPeaks -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -mask planck14/HemAsym_north_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot F -out peaks/nasym_cold_nu0p5.dat
./Stack -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -peaks peaks/sasym_cold_nu0p5.dat -field T -mask planck14/HemAsym_south_int_mask_120a_0128.fits -randrot T -out stacked/sasym_120a_cold_nu0p5
./Stack -map planck14/dx11_v2_smica_int_cmb_120a_0128_hp_10_20.fits -peaks peaks/nasym_cold_nu0p5.dat -field T -mask planck14/HemAsym_north_int_mask_120a_0128.fits -randrot T -out stacked/nasym_120a_cold_nu0p5
fi
   
for i in `seq 29 99`
do
    ./GetPeaks -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -mask planck14/HemAsym_south_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot T -out peaks/ffp8_${i}_sasym_hot_nu0p5.dat
    ./GetPeaks -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -mask planck14/HemAsym_north_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot T -out peaks/ffp8_${i}_nasym_hot_nu0p5.dat
    ./Stack -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -peaks peaks/ffp8_${i}_sasym_hot_nu0p5.dat -field T -mask planck14/HemAsym_south_int_mask_120a_0128.fits -randrot T -out stacked/ffp8_${i}_sasym_120a_hot_nu0p5
    ./Stack -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -peaks peaks/ffp8_${i}_nasym_hot_nu0p5.dat -field T -mask planck14/HemAsym_north_int_mask_120a_0128.fits -randrot T -out stacked/ffp8_${i}_nasym_120a_hot_nu0p5
    ./GetPeaks -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -mask planck14/HemAsym_south_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot F -out peaks/ffp8_${i}_sasym_cold_nu0p5.dat
    ./GetPeaks -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -mask planck14/HemAsym_north_int_mask_120a_0128.fits -peak RANDOM -orient RANDOM -nu 0.5  -hot F -out peaks/ffp8_${i}_nasym_cold_nu0p5.dat
    ./Stack -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -peaks peaks/ffp8_${i}_sasym_cold_nu0p5.dat -field T -mask planck14/HemAsym_south_int_mask_120a_0128.fits -randrot T -out stacked/ffp8_${i}_sasym_120a_cold_nu0p5
    ./Stack -map ffp8/ffp8_smica_int_${i}_120a_0128_hp_10_20.fits -peaks peaks/ffp8_${i}_nasym_cold_nu0p5.dat -field T -mask planck14/HemAsym_north_int_mask_120a_0128.fits -randrot T -out stacked/ffp8_${i}_nasym_120a_cold_nu0p5    
done    





