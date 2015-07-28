#!/bin/bash
mapname=dust
nside=0
radius=0.
if [ ${nside} == "0" ]
then
    mask=dust/lat30_mask_n1024.fits
else
    mask=dust/lat30_mask_n${nside}.fits
fi
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_15a.fits -prefix gauss_15a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_30a.fits -prefix gauss_30a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}
./MM0 -map dust/${mapname}_i_n1024_gauss_sim_60a.fits -prefix gauss_60a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}
./MM0 -map dust/${mapname}_i_n1024_15a.fits -prefix dust_15a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}
./MM0 -map dust/${mapname}_i_n1024_30a.fits -prefix dust_30a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}
./MM0 -map dust/${mapname}_i_n1024_60a.fits -prefix dust_60a -numin -2.1 -numax 2.1 -nnu 22 -nside ${nside} -mask ${mask} -radius ${radius}

#map2gif -inp original_map.fits -out original_map.gif -bar T -min -420 -max 420


