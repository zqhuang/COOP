/*
Documentation:
This asymptote script makes a general 2d plots from a text file.
The format of the text file is:
=========================================================
=========================================================
figure width and figure height in inches (2 real numbers)
Figure caption (string)
x label     (string)
y label     (string)
xlog ylog zlog (3 integers, 0 means linear scale, otherwise log scale)
clip  (integer, if 0 no clipping, if !=0 clip the figure in the range given below)
x_min x_max num_xticks(2 real numbers setting the boundary of x range, for axes, if x_min > x_max this will be ignored; and maximum number of xticks)
y_min y_max num_yticks (2 real numbers setting the boundary of y range, if y_min > y_max this will be ignored; and maximum number of yticks)
n = number of blocks (integer)
... n blocks sections.... (each section contains an indepedent block)
=========================================================
=========================================================
More about the blocks:
---------------------------------------------------------
Each block can be either DOTS, LINES, CURVE, GAUSSIAN, LABELS, ARROWS, CONTOUR, CLIP, LEGEND, EXTRA_AXIS, DENSITY, EXPAND, ERRORBARS
---------------------------------------------------------
Format of DOTS block
---------------------------------------------------------
DOTS  (string, specify that this is a DOTS block)
n = number of dots to be ploted (integer)
color and line type (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "blue_solid_1", "cyan_dotted_1.5", "black_dashdotted_0.8", ...)
dot symbol (string, if use "dot" or "DOT" then a dot is plotted, otherwise the string itself is plotted)
coordinates of dot #1 (2 real numbers)
coordinates of dot #2 (2 real numbers)
....
coordinates of dot #n (2 real numbers)
---------------------------------------------------------
Format of LINES block
---------------------------------------------------------
LINES  (string, specify that this is a LINES block)
n = number of lines to be ploted (integer)
color and line type (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "brown_longdashed_1", "purple_longdashdotted_1.5", "black_dashdotted_0.8", ...)
coordinates of the beginning and the end of line #1 (4 real numbers)
coordinates of the beginning and the end of line #2 (4 real numbers)
....
coordinates of the beginning and the end of line #n (4 real numbers)
---------------------------------------------------------
Format of CURVE block
---------------------------------------------------------
CURVE  (string, specify that this is a CURVE block)
n = number of points to be connected (integer)
legend (string)
color and line type (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "brown_longdashed_1", "purple_longdashdotted_1.5", "black_dashdotted_0.8", ...)
marker (string; NULL for no marker)
coordinates of point #1 (2 real numbers)
coordinates of point #2 (2 real numbers)
....
coordinates of point #n (2 real numbers)
---------------------------------------------------------
Format of GAUSSIAN block
---------------------------------------------------------
GAUSSIAN  (string, specify that this is a CURVE block)
n_sigma = number of sigmas (real)
legend (string)
color and line type (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "brown_longdashed_1", "purple_longdashdotted_1.5", "black_dashdotted_0.8", ...)
marker (string; NULL for no marker)
mean standard-deviation normalization-of-peak (3 real)
---------------------------------------------------------
Format of  type I CONTOUR block (with known path)
---------------------------------------------------------
CONTOUR  (string, specify that this is a CONTOUR block)
1  (integer, specify that this is a type I contour)
color and line type to fill (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "brown_longdashed_1", "purple_longdashdotted_1.5", "black_dashdotted_0.8", ... use NULL if you do not want to fill)
color and line type for border (string, color, linetype and linewidth connected with underscore, such as "red_dashed_0.5", "brown_longdashed_1", "purple_longdashdotted_1.5", "black_dashdotted_0.8", ...use NULL if you do not want to draw the border)
smooth (integer, 1 for "do smoothing" and 0 for "no smoothing")
npath (integer)
n1 = number of points to be connected for path #1 (integer)
coordinates of point #1 of path #1 (2 real numbers)
coordinates of point #2 of path #1 (2 real numbers)
....
coordinates of point #n of path #1 (2 real numbers)
n1 = number of points to be connected for path #2 (integer)
coordinates of point #1 of path #2 (2 real numbers)
coordinates of point #2 of path #2 (2 real numbers)
....
coordinates of point #n of path #2 (2 real numbers)
......
n_npath = number of points to be connected for path #npath (integer)
coordinates of point #1 of path #npath (2 real numbers)
coordinates of point #2 of path #npath (2 real numbers)
....
coordinates of point #n of path #npath (2 real numbers)
---------------------------------------------------------
Format of  type II CONTOUR block (uniform density array)
---------------------------------------------------------
CONTOUR  (string, specify that this is a CONTOUR block)
2  (integer, specify that this is a type II contour)
xmin xmax (2 real numbers)
ymin ymax (2 real numbers)
nc  (integer, specify the number of contours)
cvals (nc real numbers specify the function values where contours are drawn)
color and line type to fill for cval[0] (string)
color and line type for border for cval[0] (string)
color and line type to fill for cval[1] (string)
color and line type for border for cval[1] (string)
...
color and line type to fill for cval[nc-1] (string)
color and line type for border for cval[nc-1] (string)
smooth (integer, 1 for "do smoothing" and 0 for "no smoothing")
nx ny (two integers specify the diension of z)
data row #1 (ny real numbers)
data row #2 (ny real numbers)
....
data row #nx (ny real numbers)
---------------------------------------------------------
Format of CLIP block
---------------------------------------------------------
CLIP  (string, specify that this is a CLIP block)
n = number of points to be connected (integer)
smooth (integer, 1 for "do smoothing" and 0 for "no smoothing")
coordinates of point #1 (2 real numbers)
coordinates of point #2 (2 real numbers)
....
coordinates of point #n (2 real numbers)
---------------------------------------------------------
Format of LABELS block
---------------------------------------------------------
LABELS  (string, specify that this is a LABELS block)
n = number of labels to be written (integer)
color and line type (string, color, linetype and linewidth connected with underscore, such as "white_solid_0.5", "darkgreen_longdashed_1", "magenta_longdashdotted_1.5", "orange_dashdotted_0.8", ...)
coordinates of label #1 (2 real numbers)
label #1  (string)
coordinates of label #2 (2 real numbers)
label #2  (string)
....
coordinates of label #n (2 real numbers)
label #n  (string)
---------------------------------------------------------
Format of ARROWS block
---------------------------------------------------------
ARROWS  (string, specify that this is a ARROWS block)
n = number of arrows to be written (integer)
color and line type (string, color, linetype and linewidth connected with underscore, such as "white_solid_0.5", "darkgreen_longdashed_1", "magenta_longdashdotted_1.5", "orange_dashdotted_0.8", ...)
x1_start y1_start  x1_end y1_end
x2_start y2_start  x2_end y2_end
...
xn_start yn_start  xn_end yn_end
---------------------------------------------------------
Format of LEGEND block
---------------------------------------------------------
LEGEND  (string, specify that this is a LEGEND block)
locaion (string, can be "N", "E", "W", "S" or "NULL")
x, y (2 real numbers specify where the legend should be put; only present if location is "NULL")
cols (integer, specify how many columns)
---------------------------------------------------------
Format of a nx by ny regular grids DENSITY block 
---------------------------------------------------------
DENSITY
label of the palette     (string)
x_min x_max  (2 real numbers for the density plot)
y_min y_max  (2 real numbers for the density plot)
z_min z_max  (2 real numbers for the density plot)
0 (here 0 means regular grids)
nx ny (2 integers specifying the size of data)
row #1 of data (ny real numbers, mapping to x = x_min and y = y_min, y_min + dy, ..., y_max)
row #2 of data (ny real numbers)
...
row #nx of data (ny real numbers, mapping to x = x_max and y = y_min, ..., y_max)
---------------------------------------------------------
Format of an irregularly sampled DENSITY block 
---------------------------------------------------------
DENSITY
x_min x_max  (2 real numbers for the density plot)
y_min y_max  (2 real numbers for the density plot)
z_min z_max  (2 real numbers for the density plot)
1 (here 1 means irregular samples)
ndata (integer, number of samples)
sample #1: x1, y1, z1 (3 real numbers)
sample #2: x2, y2, z2 (3 real numbers)
...
sample #n: x_ndata, y_ndata, z_ndata (3 real numbers)
========================================================
Format of EXTRA_AXIS
========================================================
EXTRA_AXIS (string, specify that it is an EXTRAAXIS block)
location  (string, "top" for top axis, "right" for right axis)
label (string)
log  (integer 1 for logarithm and 0 for linear)
xmin, xmax (2 real numbers)
========================================================
More about the strings
---------------------------------------------------------
NULL = empty string
If you want to turn off some labels (say, caption or x label), you can use NULL. Do not use empty lines (which will be skipped).
---------------------------------------------------------
Latex symbols quoted in $ can be used in any labels. For example you can set x label to be 
$\sqrt{x^2+1}$ 
or caption to be
measured $Q_r$
---------------------------------------------------------
# is the comment string, any line started with a single # will be ignored.
A double ## is treated as a literal #. That is to say, a line
#comment 1
is ignored. 
But
##comment ##2
is not ignored and will be read in as a string "#comment #2"
---------------------------------------------------------
End of Documentation
*/
private import math;
import graph_splinetype;
import graph_settings;
import graph;
import palette;
import contour;
//=============== global variables;
settings.outformat="pdf";

string uppercase(string str){
  string [][] table = new string[26][2];
  table[0][1] = 'A';
  table[0][0] = 'a';
  table[1][1] = 'B';
  table[1][0] = 'b';
  table[2][1] = 'C';
  table[2][0] = 'c';
  table[3][1] = 'D';
  table[3][0] = 'd';
  table[4][1] = 'E';
  table[4][0] = 'e';
  table[5][1] = 'F';
  table[5][0] = 'f';
  table[6][1] = 'G';
  table[6][0] = 'g';
  table[7][1] = 'H';
  table[7][0] = 'h';
  table[8][1] = 'I';
  table[8][0] = 'i';
  table[9][1] = 'J';
  table[9][0] = 'j';
  table[10][1] = 'K';
  table[10][0] = 'k';
  table[11][1] = 'L';
  table[11][0] = 'l';
  table[12][1] = 'M';
  table[12][0] = 'm';
  table[13][1] = 'N';
  table[13][0] = 'n';
  table[14][1] = 'O';
  table[14][0] = 'o';
  table[15][1] = 'P';
  table[15][0] = 'p';
  table[16][1] = 'Q';
  table[16][0] = 'q';
  table[17][1] = 'R';
  table[17][0] = 'r';
  table[18][1] = 'S';
  table[18][0] = 's';
  table[19][1] = 'T';
  table[19][0] = 't';
  table[20][1] = 'U';
  table[20][0] = 'u';
  table[21][1] = 'V';
  table[21][0] = 'v';
  table[22][1] = 'W';
  table[22][0] = 'w';
  table[23][1] = 'X';
  table[23][0] = 'x';
  table[24][1] = 'Y';
  table[24][0] = 'y';
  table[25][1] = 'Z';
  table[25][0] = 'z';
  return replace(str, table);
}

string lowercase(string str){
  string [][] table = new string[26][2];
  table[0][0] = 'A';
  table[0][1] = 'a';
  table[1][0] = 'B';
  table[1][1] = 'b';
  table[2][0] = 'C';
  table[2][1] = 'c';
  table[3][0] = 'D';
  table[3][1] = 'd';
  table[4][0] = 'E';
  table[4][1] = 'e';
  table[5][0] = 'F';
  table[5][1] = 'f';
  table[6][0] = 'G';
  table[6][1] = 'g';
  table[7][0] = 'H';
  table[7][1] = 'h';
  table[8][0] = 'I';
  table[8][1] = 'i';
  table[9][0] = 'J';
  table[9][1] = 'j';
  table[10][0] = 'K';
  table[10][1] = 'k';
  table[11][0] = 'L';
  table[11][1] = 'l';
  table[12][0] = 'M';
  table[12][1] = 'm';
  table[13][0] = 'N';
  table[13][1] = 'n';
  table[14][0] = 'O';
  table[14][1] = 'o';
  table[15][0] = 'P';
  table[15][1] = 'p';
  table[16][0] = 'Q';
  table[16][1] = 'q';
  table[17][0] = 'R';
  table[17][1] = 'r';
  table[18][0] = 'S';
  table[18][1] = 's';
  table[19][0] = 'T';
  table[19][1] = 't';
  table[20][0] = 'U';
  table[20][1] = 'u';
  table[21][0] = 'V';
  table[21][1] = 'v';
  table[22][0] = 'W';
  table[22][1] = 'w';
  table[23][0] = 'X';
  table[23][1] = 'x';
  table[24][0] = 'Y';
  table[24][1] = 'y';
  table[25][0] = 'Z';
  table[25][1] = 'z';
  return replace(str, table);
}

string trim_string(string rawstr){
   int istart = 0;
   int iend = length(rawstr) - 1;
   while (substr(rawstr, istart, 1) == " " || substr(rawstr, istart, 1) == "\t"
 || substr(rawstr, istart, 1) == "\b" || substr(rawstr, istart, 1) == "\v" || substr(rawstr, istart, 1) == "\n" || substr(rawstr, istart, 1) == "\r") {
       ++istart;
       if(istart > iend) return "";}
   while (substr(rawstr, iend, 1) == " " || substr(rawstr, iend, 1) == "\t"
 || substr(rawstr, iend, 1) == "\b" || substr(rawstr, iend, 1) == "\v" || substr(rawstr, iend, 1) == "\n" || substr(rawstr, iend, 1) == "\r") {
       --iend;
       if(iend < istart) return "";}
   return substr(rawstr, istart, iend-istart+1);}


string nostring(real x){
  return "";}



ticks LeftTicksNoLabel(Label format="", ticklabel ticklabel=nostring,
		       bool beginlabel=true, bool endlabel=true,
		       int N=0, int n=0, real Step=0, real step=0,
		       bool begin=true, bool end=true, tickmodifier modify=None,
		       real Size=0, real size=0, bool extend=false,
		       pen pTick=nullpen, pen ptick=nullpen)
{
  return Ticks(-1,format,ticklabel,beginlabel,endlabel,N,n,Step,step,
	       begin,end,modify,Size,size,extend,pTick,ptick);
}

ticks LeftTicksNoLabel = LeftTicksNoLabel();


ticks RightTicksNoLabel(Label format = "", ticklabel ticklabel=nostring,
			bool beginlabel=true, bool endlabel=true,
			int N=0, int n=0, real Step=0, real step=0,
			bool begin=true, bool end=true, tickmodifier modify=None,
			real Size=0, real size=0, bool extend=false,
			pen pTick=nullpen, pen ptick=nullpen)
{
  return Ticks(1,format,ticklabel,beginlabel,endlabel,N,n,Step,step,
	       begin,end,modify,Size,size,extend,pTick,ptick);
}

ticks RightTicksNoLabel = RightTicksNoLabel();
 
axis YEqualsCenter(real y, bool extend=true){   
  return new void(picture pic, axisT axis){
    axis.type=Value;
    axis.value=pic.scale.y.T(y);
    axis.position=0.5;
    axis.side=right;
    axis.align=S;
    axis.extend=extend;
  };}
 
axis XEqualsCenter(real x, bool extend=true){ 
  return new void(picture pic, axisT axis){
    axis.type=Value;
    axis.value=pic.scale.x.T(x);
    axis.position=0.5;
    axis.side=left;
    axis.align=W;
    axis.extend=extend;
  };}


axis XEqualsRight(real x, bool extend=false){ 
  return new void(picture pic, axisT axis){
    axis.type=Value;
    axis.value=pic.scale.x.T(x);
    axis.position=0.5;
    axis.side=right;
    axis.align=E;
    axis.extend=extend;
  };}
 

axis YEqualsTop(real y, bool extend=false){   
  return new void(picture pic, axisT axis){
    axis.type=Value;
    axis.value=pic.scale.y.T(y);
    axis.position=0.5;
    axis.side=left;
    axis.align=N;
    axis.extend=extend;
  };}

//======= string and io functions ======



string fetch_string(file fin){
  if(eof(fin)) return "END_OF_FILE";
  string getstr = fin;
  getstr = trim_string(getstr);
  int nlines = 0;
  while(getstr == "" || substr(getstr, 0, 1) == "#" ){
    if(eof(fin)) return "END_OF_FILE";
    getstr = fin; 
    ++nlines;
    if(nlines > 100) abort("Too many comment lines?");}
  if(getstr=="NULL") return "";
  return getstr;}
 

// string__TRANSFORM__transformations
Label fetch_label(file fin){
  string fullstr = fetch_string(fin);
  string sbreak[] = split(fullstr, "__TRANSFORM__");
  if(sbreak.length == 0 || sbreak.length > 2) return Label("");
  if(sbreak.length == 1) return Label(s = fullstr);
  Label l = Label(sbreak[0]);
  string ss[] = split(sbreak[1], ":");   
  for ( int i = 0; i<ss.length; ++i){
    if(substr(ss[i], 0, 1) == "R"){  //rotation
      l = rotate( ((real) substr(ss[i],1)) ) * l; }
    else if(substr(ss[i],0,1) == "X"){ //xscale
      l = xscale( ((real) substr(ss[i],1)) ) * l;}
    else if(substr(ss[i],0,1) == "Y"){ //yscale
      l = yscale( ((real) substr(ss[i],1)) ) * l;}
    else if(substr(ss[i],0,1) == "S"){ //x and y scale
      l = scale( ((real) substr(ss[i],1)) ) * l;}	   
    else if(substr(ss[i],0,1) == "T"){ //x and y scale
      string xy[] = split(substr(ss[i],1), "|");
      l = shift( ((real) xy[0]), ((real) xy[1]))*l;} }
  return l;
}
 
pen whitepen_from_string(string fullstr){
  string sbreak[] = split(fullstr, "_");
  if(sbreak.length == 0) return currentpen;
  string cstr;
  pen colorpen;
  colorpen = white;
  if(sbreak.length == 1) return colorpen;
  cstr = trim_string(sbreak[1]);
  if(cstr == "dotted" || cstr=="dot")
    colorpen = colorpen + dotted;
  else if(cstr == "dashed" || cstr=="dash")
    colorpen = colorpen + dashed;
  else if(cstr == "longdashed" || cstr == "longdash")
    colorpen = colorpen + longdashed;
  else if(cstr == "dashdotted" || cstr == "dotdashed" || cstr == "dashdot" || cstr == "dotdash")
    colorpen = colorpen + dashdotted;
  else if(cstr == "longdashdotted" || cstr == "longdotdashed" || cstr == "longdotdash" || cstr == "longdashdot")
    colorpen = colorpen + longdashdotted;
  else
    colorpen = colorpen + solid;
  if(sbreak.length == 2) return colorpen;
  real wid = (real) sbreak[2];
  colorpen = colorpen + wid;
  return colorpen;}
 
pen pen_from_string(string fullstr){
  string sbreak[] = split(fullstr, "_");
  if(sbreak.length == 0) return currentpen;
  string cstr = lowercase(trim_string(sbreak[0]));
  pen colorpen;
  if(cstr == "red" || cstr == "r")
    colorpen = red;
  else if( cstr == "blue" || cstr == "b")
    colorpen = blue;
  else if( cstr == "green" || cstr == "g")
    colorpen = green;
  else if( cstr == "black" || cstr == "k")
    colorpen = black;
  else if( cstr == "white" || cstr == "w")
    colorpen = white;
  else if( cstr == "magenta" || cstr == "m")
    colorpen = magenta;
  else if(cstr == "violet" || cstr == "v")
    colorpen = rgb(0.55, 0.22, 0.79);
  else if(cstr == "gold")
    colorpen = rgb(0.83, 0.63, 0.09);
  else if(cstr == "maroon")
    colorpen = rgb(0.51, 0.02, 0.25);
  else if( cstr == "cyan" ||  cstr == "c" || cstr == "turquoise")
    colorpen = cyan;
  else if( cstr == "orange" || cstr == "o")
    colorpen = orange;
  else if( cstr == "purple" || cstr == "p")
    colorpen = purple;
  else if( cstr == "brown")
    colorpen = brown;
  else if( cstr == "gray" || cstr == "grey")
    colorpen = gray;
  else if( cstr == "yellow" || cstr == "y")
    colorpen = yellow;
  else if( cstr == "olive")
    colorpen = olive;
  else if( cstr == "pink")
    colorpen = pink;
  else if( cstr == "lightgray" || cstr == "lightgrey")
    colorpen = lightgray;
  else if( cstr == "darkgray" || cstr == "darkgrey")
    colorpen = darkgray;
  else if( cstr == "darkgreen")
    colorpen = darkgreen;
  else if( cstr == "darkbrown")
    colorpen = darkbrown;
  else if( cstr == "darkblue")
    colorpen = darkblue;
  else if( cstr == "darkred")
    colorpen = rgb(0.5, 0., 0.);
  else if( cstr == "darkcyan")
    colorpen = darkcyan;
  else if( cstr == "darkmagenta")
    colorpen = darkmagenta;
  else if( cstr == "springgreen")
    colorpen = springgreen;
  else if( cstr == "lawngreen")
    colorpen = rgb(0.53, 0.97, 0.10);
  else if(cstr == "skyblue")
    colorpen = rgb(0.24, 0.6, 1.);
  else if( cstr == "royalblue")
    colorpen = royalblue;
  else if( cstr == "slateblue")
    colorpen = rgb(0.21, 0.45, 0.78);
  else if(cstr == "lightred")
    colorpen = rgb(1., 0.1, 0.1);
  else if(cstr == "lightblue")
    colorpen = rgb(0.1, 0.1, 1.);
  else if(cstr == "lightgreen")
    colorpen = rgb(0.1, 1., 0.1);
  else if( cstr == "invisible")
    colorpen = invisible;
  else{
    string[] rgbstr = split(cstr,":");
    string genre = lowercase(trim_string(rgbstr[0]));
    if(rgbstr.length == 5 && (genre == 'rgba' || genre == 'RGBA')){
      real rval = ((real) rgbstr[1])/255.;
      real gval = ((real) rgbstr[2])/255.;
      real bval = ((real) rgbstr[3])/255.;
      colorpen = rgb(rval, gval, bval) + opacity( ((real) rgbstr[4]) );
    }
    else if(rgbstr.length == 4 &&  (genre == "rgb" || genre == "RGB")){
      real rval = ((real) rgbstr[1])/255.;
      real gval = ((real) rgbstr[2])/255.;
      real bval = ((real) rgbstr[3])/255.;
      colorpen = rgb(rval, gval, bval);}
    else if(rgbstr.length == 2 && (genre == "HEX" || genre == "hex")){
      colorpen = rgb(trim_string(rgbstr[1]));}
    else if(rgbstr.length == 5 && (genre == "cmyk" || genre == "CMYK")){
      real cval = ((real) rgbstr[1])/255.; 
      real mval = ((real) rgbstr[2])/255.;
      real yval = ((real) rgbstr[3])/255.;
      real kval = ((real) rgbstr[4])/255.;
      colorpen = cmyk(cval, mval, yval, kval);}
    else if(rgbstr.length == 2 && (genre == "gray" || genre == "GRAY" || genre == "grey" || genre == "GREY")){
      real grayval = ((real) rgbstr[1])/255.;
      colorpen = gray(grayval); }
    else
      colorpen = currentpen; }
  if(sbreak.length == 1) return colorpen;
  cstr = trim_string(sbreak[1]);
  if(cstr == "dotted" || cstr=="dot")
    colorpen = colorpen + dotted;
  else if(cstr == "dashed" || cstr=="dash")
    colorpen = colorpen + dashed;
  else if(cstr == "longdashed" || cstr == "longdash")
    colorpen = colorpen + longdashed;
  else if(cstr == "dashdotted" || cstr == "dotdashed" || cstr == "dashdot" || cstr == "dotdash")
    colorpen = colorpen + dashdotted;
  else if(cstr == "longdashdotted" || cstr == "longdotdashed" || cstr == "longdotdash" || cstr == "longdashdot")
    colorpen = colorpen + longdashdotted;
  else
    colorpen = colorpen + solid;
  if(sbreak.length == 2) return colorpen;
  real wid = (real) sbreak[2];
  colorpen = colorpen + wid;
  if(sbreak.length == 3) return colorpen;
  real alpha = (real) sbreak[3];
  colorpen = colorpen + opacity(alpha);
  return colorpen;}

//============= color ======
pen rgb255(real r, real g, real b){
  return rgb(r/255., g/255., b/255.);
}


pen rgbint(int r, int g, int b){
  return rgb(r/255., g/255., b/255.);
}


picture make_picture(file inputfile){

 picture mypic;

 real cxmin, cxmax, cymin, cymax, czmin, czmax, aymin, aymax, axmin, axmax, azmin, azmax;
 int opaxis[] = new int[3];
 real infty = 0.99e30;
 int  topaxis = 0;
 int  rightaxis = 0;
 int num_xticks, num_yticks;
 real topaxis_xmin, topaxis_xmax, rightaxis_ymin, rightaxis_ymax;
 bool xlog, ylog, zlog, topaxis_xlog, rightaxis_ylog, doclip;
 bool xmin_adjust, xmax_adjust, ymin_adjust, ymax_adjust, zmin_adjust, zmax_adjust;
 string xlabel, ylabel, caption, topaxis_label, rightaxis_label;
 pen coorpen = black + solid + 1.5;
//=================== axis functions ===========

//=============  coordinate ========
 real xcoor(real x){
   if(xlog) 
     return log10(x);
   else 
     return x;}

 real ycoor(real y){
   if(ylog)
     return log10(y);
   else
     return y;}

 real[] read_xminxmax(file fin){
   real t[] = fin.dimension(2);
   if(t[0]<axmin) axmin = t[0];
   if(t[1]>axmax) axmax = t[1];
   return t;}

 real[] read_yminymax(file fin){
   real t[] = fin.dimension(2);
   if(t[0]<aymin) aymin = t[0];	
   if(t[1]>aymax) aymax = t[1];
   return t;}
	 

 real[] read_xy(file fin){
   real t[] = fin.dimension(2);
   if(t[0]<axmin) axmin = t[0];
   if(t[0]>axmax) axmax = t[0];
   if(t[1]<aymin) aymin = t[1];
   if(t[1]>aymax) aymax = t[1];
   return t;}

 real[] read_xyerrors(file fin){
   real t[] = fin.dimension(6);
   real fac = 1.1;
   if(t[0]-t[3]*fac<axmin) axmin = t[0] - t[3]*fac;
   if(t[0]-t[3]*fac>axmax) axmax = t[0] - t[3]*fac;
   if(t[0]+t[2]*fac>axmax) axmax = t[0] + t[2]*fac;
   if(t[0]+t[2]*fac<axmin) axmin = t[0] + t[2]*fac;
   if(t[1]-t[5]*fac<aymin) aymin = t[1] - t[5]*fac;
   if(t[1]-t[5]*fac>aymax) aymax = t[1] - t[5]*fac;
   if(t[1]+t[4]*fac>aymax) aymax = t[1] + t[4]*fac;
   if(t[1]+t[4]*fac<aymin) aymin = t[1] + t[4]*fac;
   return t;}


 real[] read_xyxy(file fin){
   real t[] = fin.dimension(4);
   if(t[0]<axmin) axmin = t[0];
   if(t[0]>axmax) axmax = t[0];
   if(t[2]<axmin) axmin = t[2];
   if(t[2]>axmax) axmax = t[2];
   if(t[1]<aymin) aymin = t[1];
   if(t[1]>aymax) aymax = t[1];
   if(t[3]<aymin) aymin = t[3];
   if(t[3]>aymax) aymax = t[3];
   return t;}

 real[] read_xyz(file fin){
   real t[] = fin.dimension(3);
   if(t[0]<axmin) axmin = t[0];
   if(t[0]>axmax) axmax = t[0];
   if(t[1]<aymin) aymin = t[1];
   if(t[1]>aymax) aymax = t[1];
   if(t[2]<azmin) azmin = t[2];
   if(t[2]>azmax) azmax = t[2];
   return t;}

 real[] read_mean_std_amp(file fin, real nsig){
   real t[] = fin.dimension(3);
   if(t[0]-nsig*t[1]<axmin) axmin = t[0]-nsig*t[1];
   if(t[0]+nsig*t[1]>axmax) axmax = t[0]+nsig*t[1];
   if(t[2]<aymin) aymin = t[2];
   if(t[2]>aymax) aymax = t[2];
   return t;
 }
 // =============================================================================
 // plot dots
 int plot_dots(file fin){
   int ndots = fin;
   if(ndots <= 0 || ndots > 100000){
     write(stdout, "igoring " + ((string) ndots) + " dots");
     return 0;}
   string  cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   string spotlabel = fetch_string(fin);
   real[] t;
   t = new real[2];
   if(spotlabel == "dot" || spotlabel == "DOT"){
     for (int i=0; i<ndots; ++i ) {
       t = read_xy(fin);
       dot(mypic,z = ( xcoor(t[0]), ycoor(t[1]) ), p = colorpen, filltype = Fill );}}
   else{
     for (int i=0; i<ndots; ++i ) {
       t = read_xy(fin);
       label(mypic, L = spotlabel, position = ( xcoor(t[0]), ycoor(t[1]) ), align = Center, p = colorpen);}}
   return ndots;}

 // =============================================================================
 //plot lines
 int plot_lines(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000){
     write(stdout,  "Too many lines: " + ((string) nlines) + " lines");
     return 0;}
   real[] pts;
   pts = new real[4];
   string cstr;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   for (int i = 0; i< nlines; ++i){
     pts = read_xyxy(fin);
     draw(mypic, ( xcoor(pts[0]), ycoor(pts[1]) ) -- ( xcoor(pts[2]) , ycoor(pts[3]) ) , p = colorpen ); }
   return nlines; }

 // =============================================================================
 //plot curve
 int plot_curve(file fin){
   int nlines = fin;
   string clegend;
   clegend = fetch_string(fin);
   if(nlines <= 0 || nlines > 100000){
     write(stdout,  'Too many points: ' + ((string) nlines) + ' points');
     return 0;}
   real[] pts;
   pts = new real[2];
   string cstr;
   cstr = fetch_string(fin);
   string mark = trim_string(fetch_string(fin));
   pen colorpen = pen_from_string(cstr);
   path curve;
   pts = read_xy(fin);
   curve = ( xcoor(pts[0]), ycoor(pts[1]) ) ;
   for (int i = 1; i< nlines; ++i){
     pts = read_xy(fin);
     curve = curve -- ( xcoor(pts[0]), ycoor(pts[1]) ) ;}
   if(mark!='0' && mark !=''){
     if(trim_string(clegend) != ""){
       if(mark == '*')
	 draw(mypic, curve,  colorpen, Mark[6], legend = clegend);
       else if(mark =='D' || mark == 'd')
	 draw(mypic, curve,  colorpen, MarkFill[2], legend = clegend);
       else if(mark == 'T' || mark=='t')
	 draw(mypic, curve,  colorpen, MarkFill[1], legend = clegend);
       else if(mark == 'o')
	 draw(mypic, curve,  colorpen, Mark[0], legend = clegend);
       else if(mark == '+' || mark == 'x')
	 draw(mypic, curve,  colorpen, Mark[5], legend = clegend);
       else
	 draw(mypic, curve,  colorpen, MarkFill[0], legend = clegend);}
     else{  
       if(mark == '*')
	 draw(mypic, curve,  colorpen, Mark[6]);
       else if(mark =='D' || mark == 'd')
	 draw(mypic, curve,  colorpen, MarkFill[2]);
       else if(mark == 'T' || mark=='t')
	 draw(mypic, curve,  colorpen, MarkFill[1]);
       else if(mark == 'o')
	 draw(mypic, curve,  colorpen, Mark[0]);
       else if(mark == '+' || mark == 'x')
	 draw(mypic, curve,  colorpen, Mark[5]);
       else
	 draw(mypic, curve,  colorpen, MarkFill[0]);}}
   else{
     if(trim_string(clegend) != "")
       draw(mypic, curve,  colorpen, legend = clegend);
     else  
       draw(mypic, curve, colorpen);}
   return nlines; }
//======================================================================
//-----------------plot a Gaussian distribution ------------------------
int plot_gaussian(file fin){
   real nsig = fin; //how many sigmas
   string clegend;
   clegend = fetch_string(fin); //legend
   if(nsig <= 0 || nsig > 5 ){
     write(stdout,  'Too many sigmas: ' + ((string) nsig) + ' sigma');
     return 0;}
   real[] pts;
   pts = new real[3];
   string cstr;
   cstr = fetch_string(fin);
   string mark = trim_string(fetch_string(fin));
   pen colorpen = pen_from_string(cstr);
   path curve;
   pts = read_mean_std_amp(fin, nsig); //mean, standard deviation
   real mean = pts[0];
   real std = pts[1];
   real amp = pts[2];
   real dx = std/20;
   real xmax = mean + std*nsig;
   real x = mean - std*nsig;
   curve = ( xcoor(x), ycoor(amp*exp(-((x-mean)/std)^2/2.)) ) ;
   while (x < xmax){
     x +=  dx;
     curve = curve -- ( xcoor(x), ycoor(amp*exp(-((x-mean)/std)^2/2.)) ) ; }
   if(mark!='0' && mark !=''){
     if(trim_string(clegend) != ""){
       if(mark == '*')
	 draw(mypic, curve,  colorpen, Mark[6], legend = clegend);
       else if(mark =='D' || mark == 'd')
	 draw(mypic, curve,  colorpen, MarkFill[2], legend = clegend);
       else if(mark == 'T' || mark=='t')
	 draw(mypic, curve,  colorpen, MarkFill[1], legend = clegend);
       else if(mark == 'o')
	 draw(mypic, curve,  colorpen, Mark[0], legend = clegend);
       else if(mark == '+' || mark == 'x')
	 draw(mypic, curve,  colorpen, Mark[5], legend = clegend);
       else
	 draw(mypic, curve,  colorpen, MarkFill[0], legend = clegend);}
     else{  
       if(mark == '*')
	 draw(mypic, curve,  colorpen, Mark[6]);
       else if(mark =='D' || mark == 'd')
	 draw(mypic, curve,  colorpen, MarkFill[2]);
       else if(mark == 'T' || mark=='t')
	 draw(mypic, curve,  colorpen, MarkFill[1]);
       else if(mark == 'o')
	 draw(mypic, curve,  colorpen, Mark[0]);
       else if(mark == '+' || mark == 'x')
	 draw(mypic, curve,  colorpen, Mark[5]);
       else
	 draw(mypic, curve,  colorpen, MarkFill[0]);}}
   else{
     if(trim_string(clegend) != "")
       draw(mypic, curve,  colorpen, legend = clegend);
     else  
       draw(mypic, curve, colorpen);}
   return round(nsig); }

 // =============================================================================
 //plot errorbars
 //ERRORBARS
 //num_errorbars
 //color_type_width for error bars
 //center symbol
 //color_type_width for center
 //x, y, dx_minus, dx_plus, dy_minus, dy_plus
 int plot_errorbars(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000){
     write(stdout,  'Too many points: ' + ((string) nlines) + ' points');
     return 0;}
   real[] pts;
   pts = new real[6];
   pair[] centers;
   pair[] dpf;
   pair[] dmf;
   centers = new pair[nlines];
   dpf = new pair[nlines];
   dmf = new pair[nlines];
   string cstr;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real barsize = fin;
   string center_symbol = fetch_string(fin);
   cstr = fetch_string(fin);  
   pen centerpen = pen_from_string(cstr);   
   for (int i = 0; i< nlines; ++i){
     pts = read_xyerrors(fin);
     centers[i] = ( pts[0], pts[1] );
     dpf[i] = ( pts[2], pts[4] );
     dmf[i] = ( pts[3], pts[5] );}
   errorbars(pic = mypic, z = centers, dp = dpf, dm = dmf, p = colorpen, size = barsize);
   if(center_symbol == "DOT"){
     for(int i=0; i<nlines; ++i)
       dot(mypic, z = centers[i], p = centerpen, filltype = Fill );}
   else if(center_symbol != ""){
     for(int i=0; i<nlines; ++i)   
       label(mypic, L = center_symbol, position = centers[i], align = Center, p = centerpen);}  
   return nlines; }
 // =============================================================================
 //plot contour
 int plot_contour(file fin){
   int ctype = fin;
   if ( ctype == 1 ){
     string strfill, strborder;
     strfill = fetch_string(fin);
     pen colorfill = pen_from_string(strfill);
     strborder = fetch_string(fin);
     pen colorborder = pen_from_string(strborder);
     string legend = trim_string(fetch_string(fin));
     int smooth = fin;
     int npaths = fin; 
     path[] curves;
     curves = new path[npaths];
     real[] pts;
     pts = new real[2];
     int nlines, totallines;
     totallines = 0;
     for(int ipath = 0; ipath < npaths; ++ipath){
       nlines = fin;
       totallines = totallines + nlines;
       if(nlines <= 0 || nlines > 100000){
	 write(stdout,  'Too many points: ' + ((string) nlines) + ' points');
	 return 0;}
       pts = read_xy(fin);
       curves[ipath] =( xcoor(pts[0]), ycoor(pts[1]) ) ;
       if( smooth != 0){
	 for (int i = 1; i< nlines; ++i){
	   pts = read_xy(fin);
	   curves[ipath] = curves[ipath] .. ( xcoor(pts[0]), ycoor(pts[1]) ) ;
	 }  
	 curves[ipath] = curves[ipath] .. cycle ; }
       else{
	 for (int i = 1; i< nlines; ++i){
	   pts = read_xy(fin);
	   curves[ipath] = curves[ipath] -- ( xcoor(pts[0]), ycoor(pts[1]) ) ;}
	 curves[ipath] = curves[ipath] -- cycle ;
       }
     }
     if(legend == ''){
       if(trim_string(strfill) != "")  fill(mypic, curves, colorfill);
       if(trim_string(strborder) != "")  draw(mypic, curves, colorborder, legend = legend);}
     else{
       if(trim_string(strfill) != "")  fill(mypic, curves, colorfill);
       if(trim_string(strborder) != "")  draw(mypic, curves, colorborder, legend =legend);}
     return totallines;}
   else if(ctype == 2 ){
     real[] t;
     t = new real[2];
     t = read_xminxmax(fin);    
     real xmin = t[0];
     real xmax = t[1];     
     t = read_yminymax(fin);        
     real ymin = t[0];
     real ymax = t[1];
     int nc = fin; //number of contours
     real cvals[] = fin.dimension(nc); //read in n countour values
     string[] strfill, strborder;
     strfill = new string[nc];
     strborder = new string[nc];
     for ( int i = 0;  i< nc; ++i){
       strfill[i] = fetch_string(fin);
       strborder[i] = fetch_string(fin);}
     int smooth = fin;
     int ns[] = fin.dimension(2);
     int nx = ns[0];
     int ny = ns[1];
     real z[][] = fin.dimension(nx, ny);
     if(smooth != 0){
       guide ct[][] = contour(z, (xmin,ymin), (xmax, ymax), cvals, join = operator .. ); 
       for(int i=ct.length-1; i>=0; --i){
	 if(trim_string(strfill[i]) != ''){
	   pen colorfill = pen_from_string(strfill[i]);    
	   for (int j=0; j<ct[i][:].length; ++j){
	     fill(mypic, ct[i][j]..cycle, colorfill);
	   }}
	 if(trim_string(strborder[i]) != ''){
	   pen colorborder = pen_from_string(strborder[i]);   
           draw(mypic, ct[i][:], colorborder);} 
       }}
     else{
       guide ct[][] = contour(z, (xmin,ymin), (xmax, ymax), cvals, join = operator --);     
       for(int i=ct.length-1; i>=0 ; --i){
	 if(trim_string(strfill[i]) != ''){
	   pen colorfill = pen_from_string(strfill[i]);    
	   for (int j=0; j< ct[i][:].length; ++j){
	     fill(mypic, ct[i][j]--cycle, colorfill);
	   }}
	 if(trim_string(strborder[i]) != ''){
	   pen colorborder = pen_from_string(strborder[i]);   
	   draw(mypic, ct[i][:], colorborder);} 
       }}
     return nx*ny;
   }
   else	
     return 0;
 }

 // =============================================================================
 //clip contour
 int plot_clip(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000){
     write(stdout,  'Too many points: ' + ((string) nlines) + ' points');
     return 0;}
   real[] pts;
   pts = new real[2];
   path curve;
   int smooth = fin;
   pts = read_xy(fin);
   curve =( xcoor(pts[0]), ycoor(pts[1]) ) ;
   if( smooth != 0){
     for (int i = 1; i< nlines; ++i){
       pts = read_xy(fin);
       curve = curve .. ( xcoor(pts[0]), ycoor(pts[1]) ) ;}
     curve = curve .. cycle ; }
   else{
     for (int i = 1; i< nlines; ++i){
       pts = read_xy(fin);
       curve = curve -- ( xcoor(pts[0]), ycoor(pts[1]) ) ;}
     curve = curve -- cycle ;
   }
   clip(mypic, curve);
   return nlines; }

 // =============================================================================
 //plot labels
 int plot_labels(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000)  {
     write(stdout, "Too many labels: " + ((string) nlines) + " lables");
     return 0;}
   string cstr;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real [] t;
   t = new real[2];
   Label l;
   for (int i = 0; i< nlines; ++i){
     t = read_xy(fin);
     l = fetch_label(fin);
     label(mypic, L = l, position = ( xcoor(t[0]), ycoor(t[1]) ), align = Center, p = colorpen);}
   return nlines;}

 int plot_labels_left(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000)  {
     write(stdout, "Too many labels: " + ((string) nlines) + " lables");
     return 0;}
   string cstr;
   Label l;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real [] t;
   t = new real[2];
   for (int i = 0; i< nlines; ++i){
     t = read_xy(fin);
     l = fetch_label(fin);
     label(mypic, L = l, position = ( xcoor(t[0]), ycoor(t[1]) ), align = LeftSide, p = colorpen);}
   return nlines;}

 int plot_labels_right(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000)  {
     write(stdout, "Too many labels: " + ((string) nlines) + " lables");
     return 0;}
   string cstr;
   Label l;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real [] t;
   t = new real[2];
   for (int i = 0; i< nlines; ++i){
     t = read_xy(fin);
     l = fetch_label(fin);
     label(mypic, L = l, position = ( xcoor(t[0]), ycoor(t[1]) ), align = RightSide, p = colorpen);}
   return nlines;}

 int plot_arrows(file fin){
   int nlines = fin;
   if(nlines <= 0 || nlines > 100000)  {
     write(stdout, "Too many arrows: " + ((string) nlines) + " lables");
     return 0;}
   string cstr;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real [] t;
   t = new real[4];
   for (int i = 0; i< nlines; ++i){
     t = read_xyxy(fin);
     draw(mypic, (t[0], t[1]) -- (t[2], t[3]), colorpen, Arrow);}
   return nlines;
 }

 int plot_annotate(file fin){
   string cstr;
   cstr = fetch_string(fin);
   pen colorpen = pen_from_string(cstr);
   real [] t;
   t = new real[4];
   t = read_xyxy(fin);
   cstr = fetch_string(fin);
   label(mypic, L = cstr, position = ( xcoor(t[2]), ycoor(t[3]) ), align = Center, p = colorpen);
   draw(mypic, (t[2]*0.9+t[0]*0.1, t[3]*0.9+t[1]*0.1) -- (t[0]*0.9+t[2]*0.1, t[1]*0.9+t[3]*0.1), colorpen, Arrow);
   return 1;
 }

 // =============================================================================
 //plot legends
 int plot_legend(file fin){
   string cstr;
   cstr = fetch_string(fin);
   if(trim_string(cstr) !=""){
     if(trim_string(cstr) == "VIRTUAL"){
       string l = fetch_string(fin);
       cstr = fetch_string(fin);
       pen colorpen = pen_from_string(cstr) + linecap(0);
       aymax = aymax + (aymax-aymin)*0.01;
       path  g = (axmax, aymax) .. cycle;
       draw(mypic,  g = g, p = colorpen, legend=l);
       colorpen = whitepen_from_string(cstr) + linecap(0); 	
       draw(mypic, g = g, p = colorpen);
       return 1;  }
     else{
       int cols = fin;
       if(trim_string(cstr) == "N")
	 add(dest = mypic, legend(mypic,cols), point(mypic,N), 20N, UnFill); 
       else if(trim_string(cstr) == "S")
	 add(mypic, legend(mypic,cols), point(mypic,S), 20S, UnFill); 
       else if(trim_string(cstr) == "W")
	 add(dest = mypic, legend(mypic,cols), point(mypic,W), 20W, UnFill);
       else 
	 add(dest = mypic, legend(mypic,cols), point(mypic,E), 20E, UnFill);}}
   else{
     real loc[] = fin.dimension(2);
     int cols = fin;
     add(dest = mypic, src = legend(mypic,cols), position = ( xcoor(loc[0]), ycoor(loc[1]) ), filltype = UnFill);}
   return 1;
 }

 int plot_legend_nobox(file fin){
   string cstr;
   cstr = fetch_string(fin);
   if(trim_string(cstr) !=""){
     if(trim_string(cstr) == "VIRTUAL"){
       string l = fetch_string(fin);
       cstr = fetch_string(fin);
       pen colorpen = pen_from_string(cstr) + linecap(0);
       aymax = aymax + (aymax-aymin)*0.01;
       path  g =  (axmax, aymax) .. cycle;
       draw(mypic,  g = g, p = colorpen, legend=l);
       colorpen = whitepen_from_string(cstr) + linecap(0); 	
       draw(mypic, g = g, p = colorpen);
       return 1;  }
     else{	  
       int cols = fin;
       if(trim_string(cstr) == "N")
	 add(mypic, legend(mypic,perline = cols, p = invisible), point(mypic,N), 20N, UnFill); 
       else if(trim_string(cstr) == "S")
	 add(mypic, legend(mypic,perline = cols, p = invisible), point(mypic,S), 20S, UnFill); 
       else if(trim_string(cstr) == "W")
	 add(mypic, legend(mypic,perline = cols, p = invisible), point(mypic,W), 20W, UnFill);
       else 
	 add(mypic, legend(mypic,perline = cols, p = invisible), point(mypic,E), 20E, UnFill);}}
   else{
     real loc[] = fin.dimension(2);
     int cols = fin;
     add(mypic, legend(mypic,perline = cols, p = invisible), ( xcoor(loc[0]), ycoor(loc[1]) ), UnFill);}
   return 1;
 }


 int plot_legend_advance(file fin){
   string cstr;  
   cstr = fetch_string(fin);
   pen  boxpen = pen_from_string(cstr);
   real xmargin = fin;
   real ymargin = fin;
   real linelength = fin;
   real hskip = fin;
   real vskip = fin;
   int cols = fin;  
   frame leg = legend(mypic,perline = cols, xmargin = xmargin*legendmargin, ymargin = ymargin*legendmargin, linelength = linelength*legendlinelength, hskip = hskip*legendhskip, vskip = vskip*legendvskip, p = boxpen);
   cstr = fetch_string(fin);  
   if(trim_string(cstr) !=""){
     if(trim_string(cstr) == "N")
       add(mypic, leg, point(mypic,N), 20N, UnFill); 
     else if(trim_string(cstr) == "S")
       add(mypic, leg, point(mypic,S), 20S, UnFill); 
     else if(trim_string(cstr) == "W")
       add(mypic, leg, point(mypic,W), 20W, UnFill);
     else 
       add(mypic, leg, point(mypic,E), 20E, UnFill);}
   else{
     real loc[] = fin.dimension(2);
     add(mypic, leg, ( xcoor(loc[0]), ycoor(loc[1]) ), UnFill);}
   return 1;
 }
 // =============================================================================
 //plot EXTRA_AXIS
 string plot_extra_axis(file fin){
   string loc = fetch_string(fin);
   string label = fetch_string(fin);
   int islog = fin;
   real minmax[] = fin.dimension(2);
   if(trim_string(loc) == "right" || trim_string(loc) == "RIGHT"){
     rightaxis = 1;
     rightaxis_label = label;
     rightaxis_ylog = (islog!=0);
     rightaxis_ymin = minmax[0];
     rightaxis_ymax = minmax[1];
     return "right axis";}
   else  if(trim_string(loc) == "top" || trim_string(loc) == "TOP"){
     topaxis = 1;
     topaxis_label = label;
     topaxis_xlog = (islog!=0);
     topaxis_xmin = minmax[0];
     topaxis_xmax = minmax[1];
     return "top axis";}
   else{
     write(stdout, "Unknown EXTRA_AXIS type");
     return "unknown";}
 }

 void plot_rightaxis(){
   picture q;
   if(rightaxis_ylog){
     q=secondaryY(primary = mypic, new void(picture pic) {
	 if(xlog)
	   scale(pic,Log,Log);
	 else
	   scale(pic,Linear,Log);
	 ylimits(pic, rightaxis_ymin, rightaxis_ymax);
	 yaxis(pic, rightaxis_label, XEqualsRight(cxmax), coorpen, LeftTicks("", begin=false,end=false));
       });}
   else{
     q=secondaryY(primary = mypic, new void(picture pic) {
         if(xlog)
	   scale(pic,Log,Linear);
         else
	   scale(pic,Linear,Linear);
         ylimits(pic, rightaxis_ymin, rightaxis_ymax);
	 yaxis(pic, rightaxis_label, XEqualsRight(cxmax),coorpen, LeftTicks("", begin=false,end=false));
       });}

   add(dest = mypic, src = q);
 }




 void plot_topaxis(){
   picture q;
   if(topaxis_xlog){
     q=secondaryX(primary = mypic, new void(picture pic) {
	 if(ylog)
	   scale(pic,Log,Log);
	 else
	   scale(pic,Log,Linear);
	 xlimits(pic, topaxis_xmin, topaxis_xmax);
	 xaxis(pic, topaxis_label, YEqualsTop(cymax), coorpen, RightTicks("", begin=false,end=false));
       });}
   else{
     q=secondaryX(primary = mypic, new void(picture pic) {
	 if(ylog)
	   scale(pic, Linear, Log);
	 else
	   scale(pic, Linear, Linear);
	 xlimits(pic, topaxis_xmin, topaxis_xmax);
	 xaxis(pic, topaxis_label, YEqualsTop(cymax), coorpen, RightTicks("",  begin=false, end=false));
       });}
   add(dest = mypic, src = q);}


 // =============================================================================
 //plot density
 int plot_density(file fin){
   real xmin, xmax, ymin, ymax, zmin, zmax, xmincoor, xmaxcoor, ymincoor, ymaxcoor;
   real[] t;
   pen [] p;
   p = new pen[256];
   string ctbl = fetch_string(fin);
   if(ctbl == "BWRainbow")
     p = BWRainbow(256);
   else if(ctbl == "Grayscale")
     p = Grayscale(256);
   else if(ctbl == "MyRainbow")
     p = Gradient(256, darkblue,  blue, cyan, green, yellow, orange, red, darkred); 
   else if(ctbl == "Planck")
     p = Gradient(256,rgb(0,0,255),rgb(0,13,255),rgb(0,26,255),rgb(0,40,255),rgb(0,53,255),rgb(0,66,255),rgb(0,80,255),rgb(0,93,255),rgb(0,106,255),rgb(0,119,255),rgb(0,132,255),rgb(0,144,255),rgb(0,157,255),rgb(0,170,255),rgb(0,182,255),rgb(0,195,255),rgb(0,208,255),rgb(0,221,255),rgb(30,222,250),rgb(60,224,245),rgb(91,226,241),rgb(121,228,236),rgb(151,230,232),rgb(182,232,227),rgb(212,234,223),rgb(242,236,218),rgb(255,233,201),rgb(255,226,176),rgb(255,219,151),rgb(255,213,126),rgb(255,206,100),rgb(255,199,75),rgb(255,193,50),rgb(255,186,25),rgb(255,180,0),rgb(255,167,0),rgb(255,155,0),rgb(255,142,0),rgb(255,130,0),rgb(255,117,0),rgb(255,105,0),rgb(255,92,0),rgb(255,80,0),rgb(244,69,0),rgb(226,61,0),rgb(208,52,0),rgb(190,43,0),rgb(172,34,0),rgb(154,26,0),rgb(136,17,0),rgb(118,8,0),rgb(100,0,0));
   else if(ctbl == "PlanckFreq")
     p = Gradient(256,rgbint(0,0,255),rgbint(4,8,255),rgbint(8,15,255),rgbint(13,45,255),rgbint(21,108,255),rgbint(28,171,255),rgbint(45,200,255),rgbint(65,219,255),rgbint(89,235,255),rgbint(131,237,253),rgbint(174,238,251),rgbint(200,239,249),rgbint(214,240,247),rgbint(228,240,245),rgbint(234,240,230),rgbint(240,241,215),rgbint(242,241,202),rgbint(244,240,185),rgbint(245,239,168),rgbint(247,237,151),rgbint(248,235,133),rgbint(249,225,102),rgbint(249,214,66),rgbint(249,200,36),rgbint(246,180,26),rgbint(243,161,17),rgbint(233,135,10),rgbint(219,106,5),rgbint(204,77,0),rgbint(189,59,12),rgbint(174,42,25),rgbint(157,27,32),rgbint(138,15,32),rgbint(118,2,32),rgbint(118,39,69),rgbint(123,88,116),rgbint(131,131,157),rgbint(151,151,177),rgbint(171,171,196),rgbint(184,184,210),rgbint(194,194,220),rgbint(204,204,230),rgbint(214,214,234),rgbint(224,224,239),rgbint(231,231,243),rgbint(236,236,246),rgbint(241,241,249),rgbint(244,244,251),rgbint(246,246,252),rgbint(248,248,253),rgbint(250,250,254),rgbint(252,252,255));
   else
     p = Rainbow(256);

   string zlabel = fetch_string(fin);
   t = new real[2];
   t = read_xminxmax(fin); //xmin, xmax
   xmin = t[0]; 
   xmax = t[1];
   t = read_yminymax(fin); // ymin, ymax
   ymin = t[0];
   ymax = t[1];
   t = fin.dimension(2); //zmin, zmax
   zmin = t[0];
   zmax = t[1];
   int irr = fin;
   if ( irr >= 1 ) {  // irregular points
     int ndata = fin;
     real f[][] =fin.dimension(ndata, 3);
     real x[];
     real y[];
     real z[];
     for(int i=0; i<ndata; ++i){
       x[i] = f[i][0];
       y[i] = f[i][1];
       z[i] = f[i][2];
       if(x[i] < axmin) axmin = x[i]; 
       if(x[i] > axmax) axmax = x[i]; 
       if(y[i] < aymin) aymin = y[i]; 
       if(y[i] > aymax) aymax = y[i]; 
       if(z[i] < azmin) azmin = z[i]; 
       if(z[i] > azmax) azmax = z[i]; 

     }
     if(xlog){
       for(int i=0; i<ndata; ++i){
	 x[i] = log10(x[i]);}}
     if(ylog){
       for(int i=0; i<ndata; ++i){
	 y[i] = log10(y[i]);}} 
     bounds density;
     if(zmin < zmax)
       density = image(mypic, x, y, z, Range(zmin, zmax), p);
     else
       density = image(mypic, x, y, z, Automatic, p);
     if(irr == 1){
       if(xlog)
	 palette(mypic,zlabel, density, (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/20.), ymin), (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/8.), ymax), Right, p, PaletteTicks("")); 
       else
	 palette(mypic,zlabel, density, (xmax+(xmax-xmin)/20., ymin), (xmax+(xmax-xmin)/8., ymax), Right, p, PaletteTicks("")); }
     return ndata;
   }
   else{ // regular points, irr<=0
     int nxy[] = fin.dimension(2);  // nx, ny
     int nx = nxy[0];
     int ny = nxy[1];
     real[][] z ;
     z = new real[nx][ny];
     z = fin.dimension(nx, ny);
     bounds density;
     if( zmin < zmax )
       density = image(mypic, z, Range(zmin, zmax), (xmin,ymin), (xmax, ymax), p);
     else
       density = image(mypic, z, Automatic, (xmin,ymin), (xmax, ymax), p);
     if(irr == 0){
       if(xlog)
	 palette(mypic,zlabel, density, (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/20.), ymin), (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/10.), ymax), Right, p, PaletteTicks(rotate(90)*Label())); 
       else
	 palette(mypic,zlabel, density, (xmax+(xmax-xmin)/20., ymin), (xmax+(xmax-xmin)/10., ymax), Right, p, PaletteTicks(rotate(90)*Label()));}
     return nx*ny; 
   }
 }

 //=======================================================================

 void plot_expand(file fin){
   real[] t = fin.dimension(4);
   real dx = axmax - axmin;
   real dy = aymax - aymin;
   axmin = axmin - dx*t[0];
   axmax = axmax + dx*t[1];
   aymin = aymin - dy*t[2];
   aymax = aymax + dy*t[3]; }


 //=======================================================================

 bool plot_block(file fin){
   string block = fetch_string(fin);
   int nlines;
   bool plotted = true;
   if(block == "DOTS"){
     nlines = plot_dots(fin);
     write(stdout, (string) nlines + ' dots are plotted.\n');}
   else if(block == "LINES"){
     nlines = plot_lines(fin);
     write(stdout, (string) nlines + ' lines are plotted.\n');}
   else if(block == "LABELS"){
     nlines = plot_labels(fin);
     write(stdout, (string) nlines + ' labels are plotted.\n');}
   else if(block == "LEFTLABELS"){
     nlines = plot_labels_left(fin);
     write(stdout, (string) nlines + ' labels are plotted.\n');}
   else if(block == "RIGHTLABELS"){
     nlines = plot_labels_right(fin);
     write(stdout, (string) nlines + ' labels are plotted.\n');}              
   else if(block == "ARROWS"){
     nlines = plot_arrows(fin);
     write(stdout, (string) nlines + ' arrows are plotted.\n');}
   else if(block == "ANNOTATE")
     nlines = plot_annotate(fin);
   else if(block == "CURVE"){
     nlines = plot_curve(fin);
     write(stdout, 'a curve is plotted from ' + ((string) nlines ) + ' points.\n');}
   else if(block == "GAUSSIAN"){
     nlines = plot_gaussian(fin);
     write(stdout, 'a Gaussian distribution is plotted for ' + ((string) nlines ) + ' sigmas.\n');}
   else if(block == "ERRORBARS"){
     nlines = plot_errorbars(fin);
     write(stdout, ((string) nlines ) + ' errorbars are plotted.\n');}
   else if(block == "CONTOUR"){
     nlines = plot_contour(fin);
     write(stdout, 'a contour is plotted from ' + ((string) nlines ) + ' points.\n');}
   else if(block == "CLIP"){
     nlines = plot_clip(fin);
     write(stdout, 'a contour is clipped from ' + ((string) nlines ) + ' points.\n');}
   else if(block == "DENSITY"){
     nlines = plot_density(fin);
     write(stdout, (string) nlines + ' density points are plotted.\n');}
   else if(block == "EXPAND"){
     plot_expand(fin);}
   else if(block == "LEGEND"){
     nlines = plot_legend(fin);
     write(stdout, 'legends are added. \n');}
   else if(block == "LEGEND_NOBOX"){
     nlines = plot_legend_nobox(fin);
     write(stdout, 'legends are added. \n');}              
   else if(block == "LEGEND_ADVANCE"){
     nlines = plot_legend_advance(fin);
     write(stdout, 'legends are added. \n');}       
   else if(block == "EXTRA_AXIS"){
     string added =  plot_extra_axis(fin);
     write(stdout, added + ' is added. \n');}
   else
     plotted = false;
   return plotted;}


 void plot_axes(){
   //==================== set up the coordinates ============
   real xmincoor, xmaxcoor, ymincoor, ymaxcoor;
   if(xmin_adjust && cxmin > axmin) cxmin = axmin;
   if(xmax_adjust && cxmax < axmax) cxmax = axmax;
   if(ymin_adjust && cymin > aymin) cymin = aymin;
   if(ymax_adjust && cymax < aymax) cymax = aymax;
   xmincoor = xcoor(cxmin);
   xmaxcoor = xcoor(cxmax);
   ymincoor = ycoor(cymin);
   ymaxcoor = ycoor(cymax);
   if(!xmin_adjust && !xmax_adjust && !ymin_adjust && !ymax_adjust && doclip) 
     clip(mypic, (xmincoor, ymincoor) -- (xmaxcoor, ymincoor) -- (xmaxcoor, ymaxcoor) -- (xmincoor, ymaxcoor) -- cycle );
   if(caption !=  '')
     label(mypic, caption, ( xmincoor*0.5+xmaxcoor*0.5, ymaxcoor+(ymaxcoor-ymincoor)*0.06 ) );
   if(topaxis == 0){
     xaxis(mypic,xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax, p=coorpen, ticks=LeftTicks, above=true);
     if(opaxis[0]==0)
       xaxis(mypic,"", axis=YEqualsCenter(cymax, false), xmin = cxmin, xmax = cxmax, p=coorpen, ticks=RightTicksNoLabel, above=true);
     else
       xaxis(mypic,"", axis=YEqualsCenter(cymax, false), xmin = cxmin, xmax = cxmax, p=coorpen, ticks=NoTicks, above=true);    
   }       
   else{
     xaxis(mypic,xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax,  p=coorpen, ticks=LeftTicks, above=true);
     plot_topaxis();}

   if(rightaxis == 0){
     yaxis(mypic,ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, p = coorpen, ticks=RightTicks(format = rotate(90)*Label()), above=true);
     if(opaxis[1]==0)
       yaxis(mypic,"",  axis=XEqualsCenter(cxmax, false), ymin = cymin, ymax = cymax, p=coorpen, ticks=LeftTicksNoLabel(format = rotate(90)*Label()), above=true);
     else
       yaxis(mypic,"",  axis=XEqualsCenter(cxmax, false), ymin = cymin, ymax = cymax, p=coorpen, ticks=NoTicks, above=true);
   }
   else{ 
     yaxis(mypic,ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, p=coorpen, ticks=RightTicks, above=true);
     plot_rightaxis();}

 }


 void set_scales(){
   if(zlog){
     if(xlog && ylog)
       scale(mypic, Log, Log, Log);
     else{
       if(xlog)
	 scale(mypic,Log, Linear, Log);
       else if(ylog)
	 scale(mypic,Linear, Log, Log);
       else
	 scale(mypic,Linear, Linear, Log);}}
   else{
     if(xlog && ylog)
       scale(mypic,Log, Log, Linear);
     else{
       if(xlog)
	 scale(mypic,Log, Linear, Linear);
       else if(ylog)
	 scale(mypic,Linear, Log, Linear);
       else
	 scale(mypic,Linear, Linear, Linear); }}
 }

// =============================== Main Routine===============================
// =============================================================================
//read in width and height of the figure, in unit inch
 real t[] = inputfile.dimension(2); 
 size(mypic,t[0]*inch, t[1]*inch, IgnoreAspect);
 // =============================================================================
 //read in the captioin, x label, y label
 caption = fetch_string(inputfile);
 xlabel = fetch_string(inputfile);
 ylabel = fetch_string(inputfile);
 // =========================================================
 //axis options
 opaxis = inputfile.dimension(3); 
 xlog = (opaxis[0] > 0);
 ylog = (opaxis[1] > 0);
 zlog = (opaxis[2] > 0);
 set_scales();
 //===================================
 // do clipping?
 int i = inputfile;
 doclip = (i != 0);
 //==================================================
 // read in x, y limits
 real[] t;
 t = new real [3];
 t = inputfile.dimension(3); //xmin, xmax, num_xticks; num_xticks is ignored in Asymptote; in python the automatic ticks are bad so need this parameter
 cxmin = t[0];
 cxmax = t[1];
 num_xticks = (int) t[2];
 xmin_adjust = (cxmin >= infty);
 xmax_adjust = (cxmax <= -infty);
 t = inputfile.dimension(3); // ymin, ymax, num_yticks; num_yticks is ignored in Asymptote;
 cymin = t[0];
 cymax = t[1];
 num_yticks = (int) t[2];
 ymin_adjust = (cymin >= infty);
 ymax_adjust = (cymax <= -infty);
 //here you might want to upgrade?
 czmin = infty;
 czmax = -infty;
 zmin_adjust = (czmin >= infty);
 zmax_adjust = (czmax <= infty);

 axmin = infty;
 axmax = -infty;
 aymin = infty;
 aymax = -infty;
 azmin = infty;
 azmax = -infty;

 //=================================================================
 //plot the blocks 
 int nblocks = inputfile;
 if(nblocks > 0){  //plot the first n blocks
   for(int iblock = 0; iblock < nblocks; ++iblock){
     if(! plot_block(inputfile)) break;}}
 else{ //plot all
   while(plot_block(inputfile));
 }
 // plot the axes
 if(num_xticks>=0 && num_yticks>=0) plot_axes();
 return mypic;
};


picture load_picture(string filename){
  file inputfile = input(filename);
  picture mypic = make_picture(inputfile);
  close(inputfile);
  return mypic;
  };




string pic1file, pic2file;
int npics;
real pos;
file fconf =input(name = "asyplot.config", check = false);
pic1file = fconf;
pic1file = trim_string(pic1file);
if(pic1file != ""){
  pic2file = fconf;
  pic2file = trim_string(pic2file);    
  if(pic2file !=""){
    npics = 2;
    pos = fconf;}
  else
    npics = 1;}
else{
  pic1file = getstring(prompt="Enter the 2d image text file: ");
  pic1file = trim_string(pic1file);   
  npics = 1;}
picture pic1, pic2;
pic1 = load_picture(pic1file);
if(npics >1){
  pic2 = load_picture(pic2file);
  add(dest = currentpicture, src= pic2.fit(), position = ( pos*inch, 0.), align = E, above = false);   
  add(dest = currentpicture, src= pic1.fit(), position = (0, 0.), align = W, filltype=UnFill, above = true);}
else
  add(dest= currentpicture, src = pic1.fit());








