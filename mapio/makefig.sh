#!/bin/bash
rm -f *.gif
map2gif -inp inpainted_map_8.fits -out inp_0008.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_16.fits -out inp_0016.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_32.fits -out inp_0032.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_64.fits -out inp_0064.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_128.fits -out inp_0128.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_256.fits -out inp_0256.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_512.fits -out inp_0512.gif -bar T -min -420 -max 420
map2gif -inp inpainted_map_1024.fits -out inp_1024.gif -bar T -min -420 -max 420
map2gif -inp original_map.fits -out original_map.gif -bar T -min -420 -max 420


