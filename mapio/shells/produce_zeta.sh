#! /bin/bash
./MSMAP ${1} tmpmaps/smica_mask.fits MULTIPLY ${1/.fits/_masked.fits}
map2gif -inp ${1/.fits/_masked.fits} -out ${3}zeta_Mean_from_${2}.gif -bar T -ttl "10^5 zeta" -min -35. -max 35.
map2gif -inp ${1} -out ${3}zeta_Fluc_from_${2}.gif -bar T -ttl "10^5 zeta" -min -35. -max 35. -sig 2
map2gif -inp ${1} -out ${3}zeta_MeanPlusFluc_from_${2}.gif -bar T -ttl "10^5 zeta" -min -35. -max 35. -sig 3
