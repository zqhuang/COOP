#!/bin/bash
mapname=dust
nside=0
radius=0.
numin=-2.3
numax=2.3
nnu=24
gf=F
if [ ${nside} == "0" ]
then
    mask=dust/lat30_mask_n1024.fits
else
    mask=dust/lat30_mask_n${nside}.fits
fi
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_15a.fits -prefix gauss_15a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_30a.fits -prefix gauss_30a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_60a.fits -prefix gauss_60a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}
./MM0 -map dust/${mapname}_i_n1024_15a.fits -prefix dust_15a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}
./MM0 -map dust/${mapname}_i_n1024_30a.fits -prefix dust_30a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}
./MM0 -map dust/${mapname}_i_n1024_60a.fits -prefix dust_60a -numin ${numin} -numax ${numax} -nnu ${nnu} -nside ${nside} -mask ${mask} -radius ${radius} -gf ${gf}

#map2gif -inp original_map.fits -out original_map.gif -bar T -min -420 -max 420


