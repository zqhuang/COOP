#! /bin/bash
maxval =  60
minval = -60
fwhm=20
nside=0512
highpass="-hpl1 20 -hpl2 40"
outfolder=zeta256_masked
#mask=""
mask="-mask planck14/smica_mask_0${fwhm}a_${nside}.fits"
#mask="-mask planck14/mask94_${nside}.fits"
dot=1
dote=1
doupdate=1

if [ $dot -ne 0 ]
then
    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp planck14/dx11_v2_smica_int_cmb_0${fwhm}a_${nside}.fits -out ${outfolder}/t2zeta_vis.fits -fwhm ${fwhm} fluc T -action T2Z -weight vis ${mask}
    fi

    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_vis.fits -out ${outfolder}/t2zeta_vis_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}

    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp planck14/dx11_v2_smica_int_cmb_0${fwhm}a_${nside}.fits -out ${outfolder}/t2zeta_latevis.fits -fwhm ${fwhm} fluc T -action T2Z -weight latevis  ${mask}
    fi
    map2gif -inp  ${outfolder}/t2zeta_latevis.fits -out ${outfolder}/t2zeta_latevis_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_latevis.fits -out ${outfolder}/t2zeta_latevis_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_latevis.fits -out ${outfolder}/t2zeta_latevis_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}

    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp planck14/dx11_v2_smica_int_cmb_0${fwhm}a_${nside}.fits -out ${outfolder}/t2zeta_recomb_slice.fits -fwhm ${fwhm} fluc T -action T2Z -weight recomb_slice  ${mask}
    fi
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/t2zeta_recomb_slice.fits -out ${outfolder}/t2zeta_recomb_slice_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}

fi

if [ $dote -ne 0 ]
then
    if [ $doupdate -ne 0 ]
    then
	./TE2Z -inp planck14/smica_TEB_0${fwhm}a_${nside}.fits -out ${outfolder}/te2zeta_vis.fits -fwhm ${fwhm} fluc T -action TE2Z -weight vis ${highpass}  ${mask}
    fi
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_vis.fits -out ${outfolder}/te2zeta_vis_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}

    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp planck14/smica_TEB_0${fwhm}a_${nside}.fits -out ${outfolder}/te2zeta_latevis.fits -fwhm ${fwhm} fluc T -action TE2Z -weight latevis ${highpass}  ${mask}
    fi
    map2gif -inp  ${outfolder}/te2zeta_latevis.fits -out ${outfolder}/te2zeta_latevis_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_latevis.fits -out ${outfolder}/te2zeta_latevis_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_latevis.fits -out ${outfolder}/te2zeta_latevis_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}

    if [ $doupdate -ne 0 ]
    then    
	./TE2Z -inp planck14/smica_TEB_0${fwhm}a_${nside}.fits -out ${outfolder}/te2zeta_recomb_slice.fits -fwhm ${fwhm} fluc T -action TE2Z -weight recomb_slice ${highpass}  ${mask}
    fi
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean.gif -sig 1 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean_plus_fluc1.gif -sig 3 -bar T -min ${minval} -max ${maxval}
    map2gif -inp  ${outfolder}/te2zeta_recomb_slice.fits -out ${outfolder}/te2zeta_recomb_slice_mean_plus_fluc2.gif -sig 5 -bar T -min ${minval} -max ${maxval}
fi
   
cd ${outfolder}
for i in `ls *.gif`
do
    convert ${i} ${i/\.gif/\.png}
done
rm -f *.gif
cd ..
