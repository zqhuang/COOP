import sys, os, string, math, re

prefix = "ffp8"

imap = "ffp8/ffp8_smica_int_00001_010a_1024.fits"
imask = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
polmap = "ffp8/ffp8_smica_pol_case3_00001_010a_1024.fits"
polmask = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
fwhm_in = 10
unit = "K"   # ffp8 use Kelvin

zetamap = imap.replace(".fits", "_converted_to_zeta.fits")
qtutmap = imap.replace(".fits", "_converted_to_TQTUT.fits")
emap = polmap.replace(".fits", "_converted_to_EB_E.fits")
bmap = polmap.replace(".fits", "_converted_to_EB_B.fits")

if(not os.path.isfile(zetamap) or not os.path.isfile(qtutmap) or not os.path.isfile(emap) or not os.path.isfile(bmap)):
    os.system("./ByProd " + imap + " " + polmap + " " + str(fwhm_in) + " " + imask + " " + polmask )

threshold = 0
fwhm = 15

check_files = True

spots_dir = "spots/"
stack_dir = "stacked/"

def getspots(map, st):
    if(threshold < 6):
        output =  prefix + "_" + st + "_threshold" + str(threshold) + "_fwhm" + str(fwhm) + ".txt"
    else:
        output =  prefix + "_" + st + "_NoThreshold_fwhm" + str(fwhm) + ".txt"    
    if(check_files):
        if( os.path.isfile(spots_dir + output)):
            print spots_dir + output + " already exists, skipping..."
            return output
    os.system("./GetSpots " + map + " " + st + " " + imask + " " + polmask + " " + prefix + " " + str(threshold) + " " + str(fwhm) + " " + str(fwhm_in))
    return output


def stack(map, spots, st):        
    if(check_files):
        if(string.find(map, "QTUT")!=-1 and (st == "QU" or st == "QrUr")):
            if(st == "QrUr"):
                output = "QTr_on_" + spots
            else:            
                output = "QT_on_" + spots
        else:
            if(st == "QrUr"):
                output = "Qr_on_" + spots
            elif (st == "zeta"):
                output = "zeta_on_" + spots
            else:
                output = st[0] + "_on_"+spots
        if(os.path.isfile(stack_dir + output)):
            print stack_dir+output+ "  already exists, skipping..."
            return
    os.system("./Stack " + map + " " + spots_dir + spots + " " + st + " " + imask + " " + polmask + " " + unit)
    os.system("../utils/fasy.sh "+stack_dir+output)


spots_tmax = getspots(imap, "Tmax")
spots_emax = getspots(emap, "Emax")
spots_bmax = getspots(bmap, "Bmax")
spots_zetamax = getspots(zetamap, "zetamax")
spots_tmax_orient = getspots(qtutmap, "Tmax_QTUTOrient")
spots_ptmax = getspots(qtutmap, "PTmax")
spots_pmax = getspots(polmap, "Pmax")


stack(imap, spots_tmax, "T")
stack(emap, spots_tmax, "E")
stack(bmap, spots_tmax, "B")
stack(zetamap, spots_tmax, "zeta")
stack(polmap, spots_tmax, "QrUr")
stack(polmap, spots_tmax, "QU")
stack(qtutmap, spots_tmax, "QrUr")
stack(qtutmap, spots_tmax, "QU")

stack(emap, spots_emax, "E")
stack(imap, spots_emax, "T")
stack(bmap, spots_emax, "B")

stack(imap, spots_bmax, "T")
stack(emap, spots_bmax, "E")
stack(bmap, spots_bmax, "B")

stack(zetamap, spots_zetamax, "zeta")
stack(imap, spots_zetamax, "T")

stack(imap, spots_tmax_orient, "T")
stack(emap, spots_tmax_orient, "E")
stack(zetamap, spots_tmax_orient, "zeta")
stack(polmap, spots_tmax_orient, "QU")
stack(qtutmap, spots_tmax_orient, "QU")

stack(polmap, spots_ptmax, "QU")
stack(polmap, spots_pmax, "QU")
stack(imap, spots_ptmax, "T")
stack(zetamap, spots_ptmax, "zeta")
stack(emap, spots_ptmax, "E")
stack(qtutmap, spots_ptmax, "QU")







