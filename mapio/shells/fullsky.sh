#!/bin/bash
#./Stack -map lowl/commander_inp_n0016.fits -peaks peaks/planck_inp.dat -field T -mask NONE -radius 180. -res 100 -min -20 -max 20 -randrot F  -out stacked/commander_inpainted_hotT_NSalign_n0016_440a -colortable Planck
#./Stack -map lowl/commander_inp_n0016.fits -peaks peaks/cold_inp.dat -field T -mask NONE -radius 180. -res 100 -min -20 -max 20 -randrot F  -out stacked/commander_inpainted_coldT_NSalign_n0016_440a -colortable Planck
./Stack -map lowl/commander_dx11d2_extdata_temp_cmb_n0016_440arc_v1_cr.fits -peaks peaks/planck_lowres.dat -field T -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -min -20 -max 20 -randrot T  -out stacked/commander_mask94_hotT_randRot_n0016_440a -colortable Planck
#./Stack -map lowl/commander_dx11d2_extdata_temp_cmb_n0016_440arc_v1_cr.fits -peaks peaks/planck_lowres.dat -field T -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -min -20 -max 20 -randrot F  -out stacked/commander_mask94_hotT_NSalign_n0016_440a -colortable Planck
./Stack -map lowl/commander_dx11d2_extdata_temp_cmb_n0016_440arc_v1_cr.fits -peaks peaks/cold_lowres.dat -field T -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -min -20 -max 20 -randrot T  -out stacked/commander_mask94_coldT_randRot_n0016_440a -colortable Planck
#./Stack -map lowl/commander_dx11d2_extdata_temp_cmb_n0016_440arc_v1_cr.fits -peaks peaks/cold_lowres.dat -field T -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -min -20 -max 20 -randrot F  -out stacked/commander_mask94_coldT_NSalign_n0016_440a -colortable Planck
#for i in `seq 0 9`
#do
#    ./Stack -map simu/simu_i_16_440a_${i}.fits -peaks peaks/planck_lowres_${i}.dat -field T -mask lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits -radius 180. -res 100 -min -20 -max 20 -randrot T  -out stacked/sim0${i}_mask94_hotT_randRot_n0016_440a -colortable Planck
#done    
    


