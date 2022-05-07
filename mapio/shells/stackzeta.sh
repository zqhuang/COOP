prefix="te2zeta"
minset="-min 1.5"
maxset="-max 8.5"
color="-colortable Rainbow"


./GetPeaks -map zetastack/${prefix}_vis_submap001.fits -out peaks/${prefix}_vis_mean.dat -hot T -orient RANDOM -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -nu 0 -peak $\\zeta$ 
./Stack -map zetastack/${prefix}_vis_submap001.fits -peaks peaks/${prefix}_vis_mean.dat -field zeta -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -radius 2 -want_pdf T -out stacked/${prefix}_vis_mean_onMeanZetaMax $minset $maxset $color

exit

./GetPeaks -map zetastack/${prefix}_vis_submap003.fits -out peaks/${prefix}_vis_realization1.dat -hot T -orient RANDOM -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -nu 0 -peak $\\zeta$ 
./Stack -map zetastack/${prefix}_vis_submap003.fits -peaks peaks/${prefix}_vis_realization1.dat -field zeta -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -radius 2 -want_pdf T -out stacked/${prefix}_vis_realization1_onZetaMax $minset $maxset $color

./GetPeaks -map zetastack/${prefix}_vis_submap005.fits -out peaks/${prefix}_vis_realization2.dat -hot T -orient RANDOM -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -nu 0 -peak $\\zeta$
./Stack -map zetastack/${prefix}_vis_submap005.fits -peaks peaks/${prefix}_vis_realization2.dat -field zeta -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -radius 2 -want_pdf T -out stacked/${prefix}_vis_realization2_onZetaMax $minset $maxset $color

./GetPeaks -map zetastack/${prefix}_vis_submap007.fits -out peaks/${prefix}_vis_realization3.dat -hot T -orient RANDOM -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -nu 0 -peak $\\zeta$
./Stack -map zetastack/${prefix}_vis_submap007.fits -peaks peaks/${prefix}_vis_realization3.dat -field zeta -mask planck14/dx11_v2_smica_int_mask_020a_0512.fits -radius 2 -want_pdf T -out stacked/${prefix}_vis_realization3_onZetaMax $minset $maxset $color


