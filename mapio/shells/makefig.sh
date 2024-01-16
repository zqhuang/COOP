#!/bin/bash
rm -f *.gif
map2gif -inp inpainted_map_8.fits -out inp_008.gif -bar T -min -300 -max 300
map2gif -inp inpainted_map_16.fits -out inp_016.gif -bar T -min -300 -max 300
map2gif -inp inpainted_map_32.fits -out inp_032.gif -bar T -min -300 -max 300
map2gif -inp inpainted_map_64.fits -out inp_064.gif -bar T -min -300 -max 300
map2gif -inp inpainted_map_128.fits -out inp_128.gif -bar T -min -300 -max 300
map2gif -inp inpainted_map_256.fits -out inp_256.gif -bar T -min -300 -max 300
map2gif -inp original_map.fits -out original_map.gif -bar T -min -300 -max 300


