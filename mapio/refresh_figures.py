import sys, os, string, math, re, glob
for fig in glob.glob("stacked/*_fwhm15.txt"):
    if( not os.path.isfile( fig.replace(".txt", ".pdf"))):
        os.system("../utils/fasy.sh " + fig )
    print fig + " done"

