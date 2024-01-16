#! /bin/bash
maxval=40
minval=-40
fwhm=20
nside=0512
outfolder=zeta${nside}
mask="-mask planck14/smica_mask_0${fwhm}a_${nside}.fits"
dot=1
dote=1
doe=1
doupdate=1
dohighpass=0
bar="-bar F"

if [ $dohighpass -ne 0 ]
then
    highpass="-hpl1 20 -hpl2 40"
    tmap="planck14/dx11_v2_smica_int_cmb_hp_20_40_0${fwhm}a_${nside}.fits"    
    emap="planck14/dx11_v2_smica_peb_case1_cmb_hp_20_40_0${fwhm}a_${nside}.fits"
    temap="planck14/smica_TEB_hp_20_40_0${fwhm}a_${nside}.fits"
else
    higpass=""
    tmap="planck14/dx11_v2_smica_int_cmb_0${fwhm}a_${nside}.fits"    
    emap="planck14/dx11_v2_smica_peb_case1_cmb_0${fwhm}a_${nside}.fits"
    temap="planck14/smica_TEB_0${fwhm}a_${nside}.fits"    
fi



if [ $dot -ne 0 ]
then
    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp ${tmap}  -out ${outfolder}/t2zeta_vis.fits -fwhm ${fwhm} fluc T -action T2Z -weight vis ${mask} ${highpass}
    fi

    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval} 
    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval} 
    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval} 


    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp ${tmap} -out ${outfolder}/t2zeta_recomb_slice.fits -fwhm ${fwhm} fluc T -action T2Z -weight recomb_slice  ${mask} ${highpass}
    fi
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval}

fi

if [ $doe -ne 0 ]
then
    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp ${emap} -out ${outfolder}/e2zeta_vis.fits -fwhm ${fwhm} fluc T -action E2Z -weight vis ${mask} ${highpass}
    fi

    map2gif -inp  ${outfolder}/e2zeta_vis.fits -out ${outfolder}/e2zeta_vis_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/e2zeta_vis.fits -out ${outfolder}/e2zeta_vis_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/e2zeta_vis.fits -out ${outfolder}/e2zeta_vis_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval}


    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp ${emap} -out ${outfolder}/e2zeta_recomb_slice.fits -fwhm ${fwhm} fluc T -action E2Z -weight recomb_slice  ${mask} ${highpass}
    fi
    map2gif -inp  ${outfolder}/e2zeta_recomb_slice.fits -out ${outfolder}/e2zeta_recomb_slice_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/e2zeta_recomb_slice.fits -out ${outfolder}/e2zeta_recomb_slice_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/e2zeta_recomb_slice.fits -out ${outfolder}/e2zeta_recomb_slice_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval}

fi


if [ $dote -ne 0 ]
then
    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp ${temap} -out ${outfolder}/te2zeta_vis.fits -fwhm ${fwhm} fluc T -action TE2Z -weight vis ${highpass}  ${mask}
    fi
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval}

    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp ${temap} -out ${outfolder}/te2zeta_recomb_slice.fits -fwhm ${fwhm} fluc T -action TE2Z -weight recomb_slice ${highpass}  ${mask}
    fi
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean.gif -sig 1 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean_plus_fluc1.gif -sig 3 ${bar} -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean_plus_fluc2.gif -sig 5 ${bar} -min ${minval} -max ${maxval}
fi
   
cd ${outfolder}
for i in `ls *.gif`
do
    convert ${i} ${i/\.gif/\.png}
done
rm -f *.gif
cd ..
