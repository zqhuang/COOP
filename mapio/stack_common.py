postfix = "_fwhm" + str(fwhm_out) + ".fits"
fullprefix = outdir + prefix

imap = fullprefix + "_I" + postfix
polmap = fullprefix + "_QU" + postfix 
tqtutmap = fullprefix + "_TQTUT" + postfix
zetamap = fullprefix + "_zeta" + postfix
zetaqzuz = fullprefix + "_zetaqzuz" + postfix
emap  = fullprefix + "_E" + postfix
bmap = fullprefix + "_B" + postfix

if(not (os.path.isfile(imap_in)  and os.path.isfile(polmap_in) and os.path.isfile(imask) and os.path.isfile(polmask))):
    if(not (os.path.isfile(imap_in))):
        print imap_in + " is missing"
    if(not (os.path.isfile(polmap_in))):
        print polmap_in + " is missing"
    if(not (os.path.isfile(imask))):
        print imask + " is missing"
    if(not (os.path.isfile(polmask))):
        print polmask + " is missing"        
    sys.exit()
    
if(not (os.path.isfile(zetamap) and os.path.isfile(tqtutmap) and os.path.isfile(emap) and os.path.isfile(bmap) and os.path.isfile(imap) and os.path.isfile(polmap) and os.path.isfile(zetaqzuz))):
    print "Producing auxiliary maps"
    os.system("./ByProd  " + fullprefix + " " + imap_in + " " + polmap_in + " " + imask + " " + polmask + " " +  str(fwhm_in) + " " + str(fwhm_out))


def strna(s):
    return filter(lambda x: x.isdigit() or x.isalpha(), s)

def getspots(inputmap, domax, peak_name, orient_name = "NULL"):
    output =  prefix  + "_fwhm" + str(fwhm_out) + "_" + strna(peak_name)    
    if(domax):
        output = output + "max"
    else:
        output = output + "min"
    output = output +  "_ORIENT" + strna(orient_name) 
    if(threshold < 10):
        output = output + "_nu"+str(threshold).replace(r".", "pt")
    else:
        output =  output 
    if(check_files):
        if( os.path.isfile(spots_dir + output+".dat")):
            print spots_dir + output + " already exists, skipping..."
            return output
    if(domax):
        col = "./GetPeaks " + inputmap + " "  + imask + " " + polmask + " T  '" + peak_name + "' '" + orient_name + "' " + spots_dir + output + " " + str(threshold)
    else:
        col = "./GetPeaks " + inputmap + " "  + imask + " " + polmask + " F  '" + peak_name + "' '" + orient_name + "' " + spots_dir + output + " " + str(threshold)
    print col
    os.system(col)
    return output


def stack(inputmap, spots, st, zmin = 1.1e31, zmax = -1.1e31, zmin2=1.1e31, zmax2 = -1.1e31):        
    if(check_files):
        output = strna(st)+ "_on_" + spots
        if(os.path.isfile(stack_dir + output + ".txt") or os.path.isfile(stack_dir + output + "_1.txt")):
            print stack_dir+output+ "  already exists, skipping..."
            return
    col = "./Stack " + inputmap + " " + imask + " " + polmask + " " +  spots_dir + spots + ".dat" + " '" + st +  "' " + stack_dir + output + " " + str(zmin) + " " + str(zmax)+ " " + str(zmin2) + " " + str(zmax2)
    print col
    os.system(col)
    if(os.path.isfile(stack_dir + output + ".txt")):
       os.system("../utils/fasy.sh "+stack_dir+output+".txt")
    if(os.path.isfile(stack_dir + output + "_1.txt")):
       os.system("../utils/fasy.sh "+stack_dir+output+"_1.txt")
    if(os.path.isfile(stack_dir + output + "_2.txt")):
       os.system("../utils/fasy.sh "+stack_dir+output+"_2.txt")


    

spots_tmax = getspots(imap, True, "$T$")
spots_tmin = getspots(imap, False, "$T$")
spots_emax = getspots(emap, True, "$E$")
spots_bmax = getspots(bmap, True, "$B$")
spots_zetamax = getspots(zetamap, True, "$\zeta$")
spots_tmax_orient = getspots(tqtutmap, True, "$T$", "$(Q_T,U_T)$")
spots_tmin_orient = getspots(tqtutmap, False, "$T$", "$(Q_T,U_T)$")
spots_zetamax_orient = getspots(zetaqzuz, True, "$\zeta$", "$(Q_\zeta,U_\zeta)$")
spots_ptmax = getspots(tqtutmap, True, "$P_T$", "$(Q_T,U_T)$")
spots_pmax = getspots(polmap, True, "$P$", "$(Q,U)$")

#stack(imap, spots_tmax, "T")
#stack(imap, spots_tmin, "T")
#stack(emap, spots_tmax, "E")
#stack(bmap, spots_tmax, "B")
#stack(zetamap, spots_tmax, "zeta")
#stack(polmap, spots_tmax, "QrUr")
#stack(polmap, spots_tmin, "QrUr")
#stack(polmap, spots_tmax, "QU")
#stack(tqtutmap, spots_tmax, "QTUT")

#stack(emap, spots_emax, "E")
#stack(imap, spots_emax, "T")
#stack(bmap, spots_emax, "B")
#stack(polmap, spots_emax, "QU")
stack(polmap, spots_emax, "QrUr")
   
#stack(imap, spots_bmax, "T")
#stack(emap, spots_bmax, "E")
#stack(bmap, spots_bmax, "B")

#stack(zetamap, spots_zetamax, "zeta")
#stack(imap, spots_zetamax, "T")
#stack(emap, spots_zetamax, "E")
#stack(bmap, spots_zetamax, "B")

stack(imap, spots_tmax_orient, "T")
stack(imap, spots_tmin_orient, "T")
#stack(emap, spots_tmax_orient, "E")
#stack(zetamap, spots_tmax_orient, "zeta")
stack(polmap, spots_tmax_orient, "QU")
stack(polmap, spots_tmin_orient, "QU")
stack(polmap, spots_tmax_orient, "QrUr")
#stack(tqtutmap, spots_tmax_orient, "QTUT")


#stack(zetamap, spots_zetamax_orient, "zeta")


stack(polmap, spots_ptmax, "QU")
stack(polmap, spots_pmax, "QU")
#stack(imap, spots_ptmax, "T")
#stack(zetamap, spots_ptmax, "zeta")
#stack(emap, spots_ptmax, "E")
#stack(tqtutmap, spots_ptmax, "QU")







