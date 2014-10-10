import sys, os, string

prefix = "simu"

imap = "simu/simulate_iqu_n2048.fits"
imask = "commander/commander_dx11d2_mask_temp_n2048_fullres_v3.fits"
polmap = "simu/simulate_iqu_10arc_n1024.fits"
polmask ="commander/commander_polmask.fits"
qtutmap = "simu/simulate_iqu_10arc_n1024_converted_to_TQTUT.fits"
emap = "simu/simulate_iqu_10arc_n1024_converted_to_TEB_submap002.fits"
bmap = "simu/simulate_iqu_10arc_n1024_converted_to__TEB_submap003.fits"
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
    os.system("./GetSpots " + map + " " + st + " " + imask + " " + polmask + " " + prefix + " " + str(threshold) + " " + str(fwhm))
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
            else:
                output = st[0] + "_on_"+spots
        if(os.path.isfile(stack_dir + output)):
            print stack_dir+output+ "  already exists, skipping..."
            return
    os.system("./Stack " + map + " " + spots_dir + spots + " " + st + " " + imask + " " + polmask)
    os.system("../utils/fasy.sh "+stack_dir+output)


spots_tmax = getspots(imap, "Tmax")
#spots_tmax_orient = getspots(imap, "Tmax_QTUTOrient")
#spots_ptmax = getspots(imap, "PTmax")
#spots_pmax = getspots(polmap, "Pmax")

#stack(imap, spots_tmax, "T")
#stack(emap, spots_tmax, "E")
#stack(polmap, spots_tmax, "QrUr")
stack(polmap, spots_tmax, "QU")
#stack(imap, spots_tmax_orient, "T")
#stack(emap, spots_tmax_orient, "E")
#stack(polmap, spots_tmax_orient, "QU")
#stack(polmap, spots_ptmax, "QU")
#stack(polmap, spots_pmax, "QU")






