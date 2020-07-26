./LLPW
./MSMAP inpainted_map_256.fits SMOOTH 440
./UDG -inp inpainted_map_256_smoothed_fwhm440arcmin.fits -out inpbar/sinp00${1}.fits -nside 16
map2gif -inp inpbar/sinp00${1}.fits -out gifs/assumed_constraint_sim00${1}_n0016_440a.gif -bar T -min -120. -max 120.
./GetPeaks -map inpbar/sinp00${1}.fits -out peaks/sinp00${1}_cold.dat -hot F -peak RANDOM -orient RANDOM -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -nu 0.
./GetPeaks -map inpbar/sinp00${1}.fits -out peaks/sinp00${1}_hot.dat -hot T -peak RANDOM -orient RANDOM -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -nu 0.
./Stack -map inpbar/sinp00${1}.fits -peaks peaks/sinp00${1}_hot.dat -field T -out stacked/assumed_constrained_sim00${1}_hot_NSalign_n0016_440a -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -colortable Planck -randrot F -min -20. -max 20.
./Stack -map inpbar/sinp00${1}.fits -peaks peaks/sinp00${1}_cold.dat -field T -out stacked/assumed_constrained_sim00${1}_cold_NSalign_n0016_440a -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -colortable Planck -randrot F -min -20. -max 20.
./Stack -map inpbar/sinp00${1}.fits -peaks peaks/sinp00${1}_hot.dat -field T -out stacked/assumed_constrained_sim00${1}_hot_randRot_n0016_440a -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -colortable Planck -randrot T -min -20. -max 20.
./Stack -map inpbar/sinp00${1}.fits -peaks peaks/sinp00${1}_cold.dat -field T -out stacked/assumed_constrained_sim00${1}_cold_randRot_n0016_440a -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -colortable Planck -randrot T -min -20. -max 20.
