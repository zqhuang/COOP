###===============================================================
### this script finds figures in your tex file and svn add them
### syntax is
### python svnaddfigs.py YourTexFileName  FigurePath
### for example
### python svnaddfigs.py oriented_stacking.tex figures
### ======   by Zhiqi Huang (zqhuang@cita.utoronto.ca) ============
import re
import os
import sys
import subprocess

##I assume you use graphicx package. For other figure packages please change the regular expression pattern accordingly.
figure_pattern = r'^\s*\\includegraphics[^\{\}]*\{([^\{\}]*)\}.*$'

###################################
def file_match(pattern, fname):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    if re.match(r'\^.+\$', pattern) :
        flag = re.I + re.M
    else:
        flag = re.I
    return re.findall(pattern, file_content, flags = flag)



pattern = figure_pattern

if(len(sys.argv) < 3):
    print "The syntax is:"
    print "python svnaddfigs.py TexFile FigurePath"
    print "TexFile is the latex file in which the script will search for figures that have been used."
    print r'FigurePath is the default figure path that you have defined in latex \graphicspath command. If you have not done so, use ./ for FigurePath.'
    sys.exit()

    
fname = sys.argv[1]
figurepath = sys.argv[2]



def add_one_figure(x):
    if x == "":
        return False
    figure =  figurepath + "/" + x
    if(os.path.isfile(figure)):
        try:
            feedback = subprocess.check_output('svn add ' + figure)
        except:
            return False
        print feedback
        return True
    figure = x
    if(os.path.isfile(figure)):
        try:
            feedback = subprocess.check_output('svn add ' + figure)
        except:
            return False
        print feedback
        return True
    print "****** " + figure + " is missing *****"
    return False
    

    
if not os.path.isfile(fname):
    print fname
    print "cannot find this file"

    
res =  file_match(pattern, fname)
anyadd = False
for x in res:
    print  "Checking " + x
    anyadd = add_one_figure(x) or anyadd


if(anyadd):
    print "New figures will be uploaded when you do svn commit"
else:
    print "No figures has been added."
