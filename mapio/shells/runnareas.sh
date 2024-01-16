#! /bin/bash
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits -lognorm T
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits -gaussianize T
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits 
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits -sim T -lognorm T
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits -sim T -gaussianize T
#./NAREAS -map planck15/dust_i.fits -mask planck15/mask_lat30.fits -sim T 
for i in `ls Commander*.txt`
do
    ../utils/fasy.sh ${i}
done

for i in `ls CMB*.txt`
do
    ../utils/fasy.sh ${i}
done
