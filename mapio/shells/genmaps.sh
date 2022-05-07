#!/bin/bash
rdeg=10
ldeg=0.
bdeg=90.
weight=vis
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}

weight=recomb_slice
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/te2zeta_${weight}.fits -out pyscripts/te2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/t2zeta_${weight}.fits -out pyscripts/t2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_mean.txt -sig 1 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_fluc1.txt -sig 3 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
./CutPatch -inp zeta0256/e2zeta_${weight}.fits -out pyscripts/e2zeta_${weight}_fluc2.txt -sig 5 -radius ${rdeg} -b ${bdeg} -l ${ldeg}
