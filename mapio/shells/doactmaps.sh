#! /bin/bash
./MSMAP act15/act15_i.fits HIGHPASS 230 270
./MSMAP act15/act15_pol.fits HIGHPASS 230 270
./MSMAP act15/act15_i_hp_230_270.fits SMOOTH 5
./MSMAP act15/act15_pol_hp_230_270.fits SMOOTH 5
