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
x_min x_max (2 real numbers setting the boundary of x range, for axes, if x_min > x_max this will be ignored)
y_min y_max (2 real numbers setting the boundary of y range, if y_min > y_max this will be ignored)
n = number of blocks (integer)
... n blocks sections.... (each section contains an indepedent block)
=========================================================
=========================================================
More about the blocks:
---------------------------------------------------------
Each block can be either DOTS, LINES, CURVE, LABELS, CONTOUR, CLIP, LEGEND, EXTRA_AXIS, or DENSITY, EXPAND
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
smooth (integer, 1 for "do smoothing" and 0 for "no smoothing")
coordinates of point #1 (2 real numbers)
coordinates of point #2 (2 real numbers)
....
coordinates of point #n (2 real numbers)
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

##by default LABELS are aligned in the center
##you can replace LABELS with LEFTLABESL or RIGHTLABELS to define the alignment
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
real cxmin, cxmax, cymin, cymax, czmin, czmax, aymin, aymax, axmin, axmax, azmin, azmax;
real infty = 0.99e30;
int  topaxis = 0;
int  rightaxis = 0;
real topaxis_xmin, topaxis_xmax, rightaxis_ymin, rightaxis_ymax;
bool xlog, ylog, zlog, topaxis_xlog, rightaxis_ylog, doclip;
bool xmin_adjust, xmax_adjust, ymin_adjust, ymax_adjust, zmin_adjust, zmax_adjust;
string xlabel, ylabel, caption, topaxis_label, rightaxis_label, textfile;
pen coorpen = black + solid + 1.5;
//=================== axis functions ===========

string nostring(real x){
       return "";}

ticks LeftTicksNoLabel(ticklabel ticklabel=nostring,
                bool beginlabel=true, bool endlabel=true,
                int N=0, int n=0, real Step=0, real step=0,
                bool begin=true, bool end=true, tickmodifier modify=None,
                real Size=0, real size=0, bool extend=false,
                pen pTick=nullpen, pen ptick=nullpen)
{
  return Ticks(-1,"",ticklabel,beginlabel,endlabel,N,n,Step,step,
               begin,end,modify,Size,size,extend,pTick,ptick);
}

ticks LeftTicksNoLabel = LeftTicksNoLabel();


ticks RightTicksNoLabel(ticklabel ticklabel=nostring,
                 bool beginlabel=true, bool endlabel=true,
                 int N=0, int n=0, real Step=0, real step=0,
                 bool begin=true, bool end=true, tickmodifier modify=None,
                 real Size=0, real size=0, bool extend=false,
                 pen pTick=nullpen, pen ptick=nullpen)
{
  return Ticks(1,"",ticklabel,beginlabel,endlabel,N,n,Step,step,
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
  string cstr = trim_string(sbreak[0]);
  pen colorpen;
  if(cstr == "red" || cstr == "RED" || cstr == "r" || cstr == "R")
      colorpen = red;
  else if( cstr == "blue" || cstr == "BLUE" || cstr == "b" || cstr == "B")
      colorpen = blue;
  else if( cstr == "green" || cstr == "GREEN" || cstr == "g" || cstr == "G")
      colorpen = green;
  else if( cstr == "black" || cstr == "BLACK")
      colorpen = black;
  else if( cstr == "white" || cstr == "WHITE" || cstr == "w" || cstr == "W")
       colorpen = white;
  else if( cstr == "magenta" || cstr == "MAGENTA" ||  cstr == "m" || cstr == "M")
      colorpen = magenta;
  else if(cstr == "violet" || cstr == "VIOLET" || cstr == "v" || cstr == "V")
      colorpen = rgb(0.55, 0.22, 0.79);
  else if(cstr == "gold" || cstr == "GOLD")
      colorpen = rgb(0.83, 0.63, 0.09);
  else if(cstr == "maroon" || cstr == "MAROON")
      colorpen = rgb(0.51, 0.02, 0.25);
  else if( cstr == "cyan" || cstr == "CYAN" ||  cstr == "turquoise" || cstr == "TURQUOISE" || cstr == "c" || cstr == "C")
      colorpen = cyan;
  else if( cstr == "orange" || cstr == "ORANGE" ||  cstr == "o" || cstr == "O")
      colorpen = orange;
  else if( cstr == "purple" || cstr == "PURPLE" || cstr == "p" || cstr == "P")
      colorpen = purple;
  else if( cstr == "brown" || cstr == "BROWN")
      colorpen = brown;
  else if( cstr == "gray" || cstr == "GRAY" || cstr == "grey" || cstr == "GREY")
      colorpen = gray;
  else if( cstr == "yellow" || cstr == "YELLOW" || cstr == "y" || cstr == "Y")
      colorpen = yellow;
  else if( cstr == "olive" || cstr == "OLIVE")
      colorpen = olive;
  else if( cstr == "pink" || cstr == "PINK")
      colorpen = pink;
  else if( cstr == "lightgray" || cstr == "LIGHTGRAY" ||  cstr == "lightgrey" || cstr == "LIGHTGREY")
      colorpen = lightgray;
  else if( cstr == "darkgray" || cstr == "DARKGRAY" || cstr == "darkgrey" || cstr == "DARKGREY")
      colorpen = darkgray;
  else if( cstr == "darkgreen" || cstr == "DARKGREEN")
      colorpen = darkgreen;
  else if( cstr == "darkbrown" || cstr == "DARKBROWN")
      colorpen = darkbrown;
  else if( cstr == "darkblue" || cstr == "DARKBLUE")
      colorpen = darkblue;
  else if( cstr == "darkred" || cstr == "DARKRED")
      colorpen = rgb(0.5, 0., 0.);
  else if( cstr == "darkcyan" || cstr == "DARKCYAN")
      colorpen = darkcyan;
  else if( cstr == "darkmagenta" || cstr == "DARKMAGENTA")
      colorpen = darkmagenta;
  else if( cstr == "springgreen" || cstr == "SPRINGGREEN")
      colorpen = springgreen;
  else if( cstr == "lawngreen" || cstr == "LAWNGREEN" || cstr == "grassgreen" || cstr == "GRASSGREEN")
      colorpen = rgb(0.53, 0.97, 0.10);
  else if(cstr == "skyblue" || cstr == "SKYBLUE")
      colorpen = rgb(0.24, 0.6, 1.);
  else if( cstr == "royalblue" || cstr == "ROYALBLUE")
      colorpen = royalblue;
  else if( cstr == "slateblue" || cstr == "SLATEBLUE")
      colorpen = rgb(0.21, 0.45, 0.78);
  else if(cstr == "lightred" || cstr == "LIGHTRED")
      colorpen = rgb(1., 0.1, 0.1);
  else if(cstr == "lightblue" || cstr == "LIGHTBLUE")
      colorpen = rgb(0.1, 0.1, 1.);
  else if(cstr == "lightgreen" || cstr == "LIGHTGREEN")
      colorpen = rgb(0.1, 1., 0.1);
  else if( cstr == "invisible" || cstr == "INVISIBLE")
      colorpen = invisible;
  else{
      string[] rgbstr = split(cstr,":");
      string gengre = trim_string(rgbstr[0]);
      if(rgbstr.length == 4 &&  (gengre == "rgb" || gengre == "RGB")){
          real rval = ((real) rgbstr[1])/255.;
	  real gval = ((real) rgbstr[2])/255.;
	  real bval = ((real) rgbstr[3])/255.;
	  colorpen = rgb(rval, gval, bval);}
      else if(rgbstr.length == 2 && (gengre == "HEX" || gengre == "hex")){
           colorpen = rgb(trim_string(rgbstr[1]));}
      else if(rgbstr.length == 5 && (gengre == "cmyk" || gengre == "CMYK")){
      	   real cval = ((real) rgbstr[1])/255.; 
           real mval = ((real) rgbstr[2])/255.;
           real yval = ((real) rgbstr[3])/255.;
           real kval = ((real) rgbstr[4])/255.;
	   colorpen = cmyk(cval, mval, yval, kval);}
      else if(rgbstr.length == 2 && (gengre == "gray" || gengre == "GRAY" || gengre == "grey" || gengre == "GREY")){
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
 return colorpen;}

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
//============= color ======
pen rgb255(real r, real g, real b){
    return rgb(r/255., g/255., b/255.);
}


pen rgbint(int r, int g, int b){
    return rgb(r/255., g/255., b/255.);
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
       dot(z = ( xcoor(t[0]), ycoor(t[1]) ), p = colorpen, filltype = Fill );}}
  else{
     for (int i=0; i<ndots; ++i ) {
       t = read_xy(fin);
       label( L = spotlabel, position = ( xcoor(t[0]), ycoor(t[1]) ), align = Center, p = colorpen);}}
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
    draw( ( xcoor(pts[0]), ycoor(pts[1]) ) -- ( xcoor(pts[2]) , ycoor(pts[3]) ) , p = colorpen ); }
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
  int smooth = fin;
  pen colorpen = pen_from_string(cstr);
  path curve;
  pts = read_xy(fin);
  curve = ( xcoor(pts[0]), ycoor(pts[1]) ) ;
  if(smooth != 0 ){
    for (int i = 1; i< nlines; ++i){
       pts = read_xy(fin);
       curve = curve .. ( xcoor(pts[0]), ycoor(pts[1]) ) ;}}
   else{
     for (int i = 1; i< nlines; ++i){
       pts = read_xy(fin);
       curve = curve -- ( xcoor(pts[0]), ycoor(pts[1]) ) ;}}
   if(trim_string(clegend) != "")
      draw(curve,  colorpen, legend = clegend);
   else  
      draw(curve, colorpen);
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
    if(trim_string(strfill) != "")  fill(curves, colorfill);
    if(trim_string(strborder) != "")  draw(curves, colorborder); 
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
        guide ct[][] = contour(z, (xmin,ymin), (xmax, ymax), cvals, join = operator .. , subsample = 1); 
       for(int i=ct.length-1; i>=0; --i){
        if(trim_string(strfill[i]) != ''){
         pen colorfill = pen_from_string(strfill[i]);    
         for (int j=0; j<ct[i][:].length; ++j){
          fill(ct[i][j]..cycle, colorfill);
	 }}
        if(trim_string(strborder[i]) != ''){
         pen colorborder = pen_from_string(strborder[i]);   
           draw(ct[i][:], colorborder);} 
    }}
    else{
       guide ct[][] = contour(z, (xmin,ymin), (xmax, ymax), cvals, join = operator -- , subsample = 1);     
      for(int i=ct.length-1; i>=0 ; --i){
       if(trim_string(strfill[i]) != ''){
         pen colorfill = pen_from_string(strfill[i]);    
         for (int j=0; j< ct[i][:].length; ++j){
           fill(ct[i][j]--cycle, colorfill);
       }}
       if(trim_string(strborder[i]) != ''){
         pen colorborder = pen_from_string(strborder[i]);   
          draw(ct[i][:], colorborder);} 
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
   clip(curve);
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
 for (int i = 0; i< nlines; ++i){
   t = read_xy(fin);
   cstr = fetch_string(fin);
   label( L = cstr, position = ( xcoor(t[0]), ycoor(t[1]) ), align = Center, p = colorpen);}
 return nlines;}

int plot_labels_left(file fin){
 int nlines = fin;
 if(nlines <= 0 || nlines > 100000)  {
    write(stdout, "Too many labels: " + ((string) nlines) + " lables");
    return 0;}
 string cstr;
 cstr = fetch_string(fin);
 pen colorpen = pen_from_string(cstr);
 real [] t;
 t = new real[2];
 for (int i = 0; i< nlines; ++i){
   t = read_xy(fin);
   cstr = fetch_string(fin);
   label( L = cstr, position = ( xcoor(t[0]), ycoor(t[1]) ), align = LeftSide, p = colorpen);}
 return nlines;}

int plot_labels_right(file fin){
 int nlines = fin;
 if(nlines <= 0 || nlines > 100000)  {
    write(stdout, "Too many labels: " + ((string) nlines) + " lables");
    return 0;}
 string cstr;
 cstr = fetch_string(fin);
 pen colorpen = pen_from_string(cstr);
 real [] t;
 t = new real[2];
 for (int i = 0; i< nlines; ++i){
   t = read_xy(fin);
   cstr = fetch_string(fin);
   label( L = cstr, position = ( xcoor(t[0]), ycoor(t[1]) ), align = RightSide, p = colorpen);}
 return nlines;}


// =============================================================================
//plot labels
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
        draw(g = g, p = colorpen, legend=l);
        colorpen = whitepen_from_string(cstr) + linecap(0); 	
	draw(g = g, p = colorpen);
	return 1;  }
     else{
       int cols = fin;
       if(trim_string(cstr) == "N")
          add(legend(cols), point(N), 20N, UnFill); 
       else if(trim_string(cstr) == "S")
          add(legend(cols), point(S), 20S, UnFill); 
       else if(trim_string(cstr) == "W")
          add(legend(cols), point(W), 20W, UnFill);
       else 
          add(legend(cols), point(E), 20E, UnFill);}}
  else{
     real loc[] = fin.dimension(2);
     int cols = fin;
     add(legend(cols), ( xcoor(loc[0]), ycoor(loc[1]) ), UnFill);}
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
    q=secondaryY(new void(picture pic) {
        if(xlog)
             scale(pic,Log,Log);
        else
             scale(pic,Linear,Log);
        ylimits(pic, rightaxis_ymin, rightaxis_ymax);
        yaxis(pic, rightaxis_label, XEqualsRight(cxmax), coorpen, LeftTicks("", begin=false,end=false));
        });}
  else{
    q=secondaryY(new void(picture pic) {
         if(xlog)
            scale(pic,Log,Linear);
         else
            scale(pic,Linear,Linear);
         ylimits(pic, rightaxis_ymin, rightaxis_ymax);
         yaxis(pic, rightaxis_label, XEqualsRight(cxmax),coorpen, LeftTicks("", begin=false,end=false));
        });}

   add(q);
}

void plot_topaxis(){
  picture q;
 if(topaxis_xlog){
    q=secondaryX(new void(picture pic) {
    if(ylog)
       scale(pic,Log,Log);
    else
       scale(pic,Log,Linear);
    xlimits(pic, topaxis_xmin, topaxis_xmax);
    xaxis(pic, topaxis_label, YEqualsTop(cymax), coorpen, RightTicks("", begin=false,end=false));
        });}
 else{
   q=secondaryX(new void(picture pic) {
   if(ylog)
       scale(pic, Linear, Log);
   else
      scale(pic, Linear, Linear);
   xlimits(pic, topaxis_xmin, topaxis_xmax);
   xaxis(pic, topaxis_label, YEqualsTop(cymax), coorpen, RightTicks("", begin=false, end=false));
        });}
 add(q);}


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
    p = Gradient(256, rgbint(0, 0, 255), rgbint(0, 2, 255), rgbint(0, 5, 255), rgbint(0, 8, 255), rgbint(0, 10, 255), rgbint(0, 13, 255), rgbint(0, 16, 255), rgbint(0, 18, 255), rgbint(0, 21, 255), rgbint(0, 24, 255), rgbint(0, 26, 255), rgbint(0, 29, 255), rgbint(0, 32, 255), rgbint(0, 34, 255), rgbint(0, 37, 255), rgbint(0, 40, 255), rgbint(0, 42, 255), rgbint(0, 45, 255), rgbint(0, 48, 255), rgbint(0, 50, 255), rgbint(0, 53, 255), rgbint(0, 56, 255), rgbint(0, 58, 255), rgbint(0, 61, 255), rgbint(0, 64, 255), rgbint(0, 66, 255), rgbint(0, 69, 255), rgbint(0, 72, 255), rgbint(0, 74, 255), rgbint(0, 77, 255), rgbint(0, 80, 255), rgbint(0, 82, 255), rgbint(0, 85, 255), rgbint(0, 88, 255), rgbint(0, 90, 255), rgbint(0, 93, 255), rgbint(0, 96, 255), rgbint(0, 98, 255), rgbint(0, 101, 255), rgbint(0, 104, 255), rgbint(0, 106, 255), rgbint(0, 109, 255), rgbint(0, 112, 255), rgbint(0, 114, 255), rgbint(0, 117, 255), rgbint(0, 119, 255), rgbint(0, 122, 255), rgbint(0, 124, 255), rgbint(0, 127, 255), rgbint(0, 129, 255), rgbint(0, 132, 255), rgbint(0, 134, 255), rgbint(0, 137, 255), rgbint(0, 139, 255), rgbint(0, 142, 255), rgbint(0, 144, 255), rgbint(0, 147, 255), rgbint(0, 150, 255), rgbint(0, 152, 255), rgbint(0, 155, 255), rgbint(0, 157, 255), rgbint(0, 160, 255), rgbint(0, 162, 255), rgbint(0, 165, 255), rgbint(0, 167, 255), rgbint(0, 170, 255), rgbint(0, 172, 255), rgbint(0, 175, 255), rgbint(0, 177, 255), rgbint(0, 180, 255), rgbint(0, 182, 255), rgbint(0, 185, 255), rgbint(0, 188, 255), rgbint(0, 190, 255), rgbint(0, 193, 255), rgbint(0, 195, 255), rgbint(0, 198, 255), rgbint(0, 200, 255), rgbint(0, 203, 255), rgbint(0, 205, 255), rgbint(0, 208, 255), rgbint(0, 210, 255), rgbint(0, 213, 255), rgbint(0, 215, 255), rgbint(0, 218, 255), rgbint(0, 221, 255), rgbint(6, 221, 254), rgbint(12, 221, 253), rgbint(18, 222, 252), rgbint(24, 222, 251), rgbint(30, 222, 250), rgbint(36, 223, 249), rgbint(42, 223, 248), rgbint(48, 224, 247), rgbint(54, 224, 246), rgbint(60, 224, 245), rgbint(66, 225, 245), rgbint(72, 225, 244), rgbint(78, 225, 243), rgbint(85, 226, 242), rgbint(91, 226, 241), rgbint(97, 227, 240), rgbint(103, 227, 239), rgbint(109, 227, 238), rgbint(115, 228, 237), rgbint(121, 228, 236), rgbint(127, 229, 236), rgbint(133, 229, 235), rgbint(139, 229, 234), rgbint(145, 230, 233), rgbint(151, 230, 232), rgbint(157, 230, 231), rgbint(163, 231, 230), rgbint(170, 231, 229), rgbint(176, 232, 228), rgbint(182, 232, 227), rgbint(188, 232, 226), rgbint(194, 233, 226), rgbint(200, 233, 225), rgbint(206, 233, 224), rgbint(212, 234, 223), rgbint(218, 234, 222), rgbint(224, 235, 221), rgbint(230, 235, 220), rgbint(236, 235, 219), rgbint(242, 236, 218), rgbint(248, 236, 217), rgbint(255, 237, 217), rgbint(255, 235, 211), rgbint(255, 234, 206), rgbint(255, 233, 201), rgbint(255, 231, 196), rgbint(255, 230, 191), rgbint(255, 229, 186), rgbint(255, 227, 181), rgbint(255, 226, 176), rgbint(255, 225, 171), rgbint(255, 223, 166), rgbint(255, 222, 161), rgbint(255, 221, 156), rgbint(255, 219, 151), rgbint(255, 218, 146), rgbint(255, 217, 141), rgbint(255, 215, 136), rgbint(255, 214, 131), rgbint(255, 213, 126), rgbint(255, 211, 121), rgbint(255, 210, 116), rgbint(255, 209, 111), rgbint(255, 207, 105), rgbint(255, 206, 100), rgbint(255, 205, 95), rgbint(255, 203, 90), rgbint(255, 202, 85), rgbint(255, 201, 80), rgbint(255, 199, 75), rgbint(255, 198, 70), rgbint(255, 197, 65), rgbint(255, 195, 60), rgbint(255, 194, 55), rgbint(255, 193, 50), rgbint(255, 191, 45), rgbint(255, 190, 40), rgbint(255, 189, 35), rgbint(255, 187, 30), rgbint(255, 186, 25), rgbint(255, 185, 20), rgbint(255, 183, 15), rgbint(255, 182, 10), rgbint(255, 181, 5), rgbint(255, 180, 0), rgbint(255, 177, 0), rgbint(255, 175, 0), rgbint(255, 172, 0), rgbint(255, 170, 0), rgbint(255, 167, 0), rgbint(255, 165, 0), rgbint(255, 162, 0), rgbint(255, 160, 0), rgbint(255, 157, 0), rgbint(255, 155, 0), rgbint(255, 152, 0), rgbint(255, 150, 0), rgbint(255, 147, 0), rgbint(255, 145, 0), rgbint(255, 142, 0), rgbint(255, 140, 0), rgbint(255, 137, 0), rgbint(255, 135, 0), rgbint(255, 132, 0), rgbint(255, 130, 0), rgbint(255, 127, 0), rgbint(255, 125, 0), rgbint(255, 122, 0), rgbint(255, 120, 0), rgbint(255, 117, 0), rgbint(255, 115, 0), rgbint(255, 112, 0), rgbint(255, 110, 0), rgbint(255, 107, 0), rgbint(255, 105, 0), rgbint(255, 102, 0), rgbint(255, 100, 0), rgbint(255, 97, 0), rgbint(255, 95, 0), rgbint(255, 92, 0), rgbint(255, 90, 0), rgbint(255, 87, 0), rgbint(255, 85, 0), rgbint(255, 82, 0), rgbint(255, 80, 0), rgbint(255, 77, 0), rgbint(255, 75, 0), rgbint(251, 73, 0), rgbint(247, 71, 0), rgbint(244, 69, 0), rgbint(240, 68, 0), rgbint(236, 66, 0), rgbint(233, 64, 0), rgbint(229, 62, 0), rgbint(226, 61, 0), rgbint(222, 59, 0), rgbint(218, 57, 0), rgbint(215, 55, 0), rgbint(211, 54, 0), rgbint(208, 52, 0), rgbint(204, 50, 0), rgbint(200, 48, 0), rgbint(197, 47, 0), rgbint(193, 45, 0), rgbint(190, 43, 0), rgbint(186, 41, 0), rgbint(182, 40, 0), rgbint(179, 38, 0), rgbint(175, 36, 0), rgbint(172, 34, 0), rgbint(168, 33, 0), rgbint(164, 31, 0), rgbint(161, 29, 0), rgbint(157, 27, 0), rgbint(154, 26, 0), rgbint(150, 24, 0), rgbint(146, 22, 0), rgbint(143, 20, 0), rgbint(139, 19, 0), rgbint(136, 17, 0), rgbint(132, 15, 0), rgbint(128, 13, 0), rgbint(125, 12, 0), rgbint(121, 10, 0), rgbint(118, 8, 0), rgbint(114, 6, 0), rgbint(110, 5, 0), rgbint(107, 3, 0), rgbint(103, 1, 0), rgbint(100, 0, 0));
  else if(ctbl == "PlanckFreq")
p = Gradient(256, rgb255(0.0, 0.00, 255.0), rgb255(0.8, 1.54, 255.0), rgb255(1.5, 3.08, 255.0), rgb255(2.3, 4.62, 255.0), rgb255(3.1, 6.15, 255.0), rgb255(3.8, 7.69, 255.0), rgb255(4.6, 9.23, 255.0), rgb255(5.4, 10.77, 255.0), rgb255(6.2, 12.31, 255.0), rgb255(6.9, 13.85, 255.0), rgb255(7.7, 15.38, 255.0), rgb255(8.5, 16.92, 255.0), rgb255(9.2, 18.46, 255.0), rgb255(10.0, 20.00, 255.0), rgb255(11.5, 32.62, 255.0), rgb255(13.1, 45.23, 255.0), rgb255(14.6, 57.85, 255.0), rgb255(16.2, 70.46, 255.0), rgb255(17.7, 83.08, 255.0), rgb255(19.2, 95.69, 255.0), rgb255(20.8, 108.31, 255.0), rgb255(22.3, 120.92, 255.0), rgb255(23.8, 133.54, 255.0), rgb255(25.4, 146.15, 255.0), rgb255(26.9, 158.77, 255.0), rgb255(28.5, 171.38, 255.0), rgb255(30.0, 184.00, 255.0), rgb255(33.8, 187.92, 255.0), rgb255(37.7, 191.85, 255.0), rgb255(41.5, 195.77, 255.0), rgb255(45.4, 199.69, 255.0), rgb255(49.2, 203.62, 255.0), rgb255(53.1, 207.54, 255.0), rgb255(56.9, 211.46, 255.0), rgb255(60.8, 215.38, 255.0), rgb255(64.6, 219.31, 255.0), rgb255(68.5, 223.23, 255.0), rgb255(72.3, 227.15, 255.0), rgb255(76.2, 231.08, 255.0), rgb255(80.0, 235.00, 255.0), rgb255(88.5, 235.31, 254.6), rgb255(97.1, 235.62, 254.2), rgb255(105.6, 235.92, 253.8), rgb255(114.2, 236.23, 253.5), rgb255(122.7, 236.54, 253.1), rgb255(131.2, 236.85, 252.7), rgb255(139.8, 237.15, 252.3), rgb255(148.3, 237.46, 251.9), rgb255(156.8, 237.77, 251.5), rgb255(165.4, 238.08, 251.2), rgb255(173.9, 238.38, 250.8), rgb255(182.5, 238.69, 250.4), rgb255(191.0, 239.00, 250.0), rgb255(193.8, 239.08, 249.6), rgb255(196.7, 239.15, 249.2), rgb255(199.5, 239.23, 248.8), rgb255(202.4, 239.31, 248.5), rgb255(205.2, 239.38, 248.1), rgb255(208.1, 239.46, 247.7), rgb255(210.9, 239.54, 247.3), rgb255(213.8, 239.62, 246.9), rgb255(216.6, 239.69, 246.5), rgb255(219.5, 239.77, 246.2), rgb255(222.3, 239.85, 245.8), rgb255(225.2, 239.92, 245.4), rgb255(228.0, 240.00, 245.0), rgb255(229.2, 240.09, 242.0), rgb255(230.4, 240.18, 239.0), rgb255(231.5, 240.27, 236.0), rgb255(232.7, 240.36, 233.0), rgb255(233.9, 240.46, 230.0), rgb255(235.1, 240.54, 227.0), rgb255(236.3, 240.64, 224.0), rgb255(237.5, 240.73, 221.0), rgb255(238.6, 240.82, 218.0), rgb255(239.8, 240.91, 215.0), rgb255(241.0, 241.00, 212.0), rgb255(241.0, 241.00, 212.0), rgb255(241.4, 240.91, 208.6), rgb255(241.7, 240.82, 205.3), rgb255(242.1, 240.73, 201.9), rgb255(242.5, 240.64, 198.5), rgb255(242.8, 240.54, 195.2), rgb255(243.2, 240.46, 191.8), rgb255(243.5, 240.36, 188.5), rgb255(243.9, 240.27, 185.1), rgb255(244.3, 240.18, 181.7), rgb255(244.6, 240.09, 178.4), rgb255(245.0, 240.00, 175.0), rgb255(245.2, 239.62, 171.5), rgb255(245.5, 239.23, 168.1), rgb255(245.7, 238.85, 164.6), rgb255(245.9, 238.46, 161.2), rgb255(246.2, 238.08, 157.7), rgb255(246.4, 237.69, 154.2), rgb255(246.6, 237.31, 150.8), rgb255(246.8, 236.92, 147.3), rgb255(247.1, 236.54, 143.8), rgb255(247.3, 236.15, 140.4), rgb255(247.5, 235.77, 136.9), rgb255(247.8, 235.38, 133.5), rgb255(248.0, 235.00, 130.0), rgb255(248.1, 232.62, 122.9), rgb255(248.3, 230.23, 115.9), rgb255(248.4, 227.85, 108.8), rgb255(248.6, 225.46, 101.8), rgb255(248.7, 223.08, 94.7), rgb255(248.9, 220.69, 87.7), rgb255(249.0, 218.31, 80.6), rgb255(249.2, 215.92, 73.5), rgb255(249.3, 213.54, 66.5), rgb255(249.5, 211.15, 59.4), rgb255(249.6, 208.77, 52.4), rgb255(249.8, 206.38, 45.3), rgb255(249.9, 204.00, 38.2), rgb255(249.3, 200.08, 36.3), rgb255(248.7, 196.15, 34.3), rgb255(248.1, 192.23, 32.4), rgb255(247.5, 188.31, 30.4), rgb255(247.0, 184.38, 28.4), rgb255(246.4, 180.46, 26.5), rgb255(245.8, 176.54, 24.5), rgb255(245.2, 172.62, 22.6), rgb255(244.6, 168.69, 20.6), rgb255(244.0, 164.77, 18.6), rgb255(243.4, 160.85, 16.7), rgb255(242.8, 156.92, 14.7), rgb255(242.2, 153.00, 12.8), rgb255(239.3, 147.12, 11.8), rgb255(236.4, 141.23, 10.8), rgb255(233.4, 135.35, 9.8), rgb255(230.5, 129.46, 8.8), rgb255(227.5, 123.58, 7.8), rgb255(224.6, 117.69, 6.9), rgb255(221.7, 111.81, 5.9), rgb255(218.7, 105.92, 4.9), rgb255(215.8, 100.04, 3.9), rgb255(212.8, 94.15, 2.9), rgb255(209.9, 88.27, 2.0), rgb255(206.9, 82.38, 1.0), rgb255(204.0, 76.50, 0.0), rgb255(201.0, 73.08, 2.5), rgb255(198.0, 69.65, 4.9), rgb255(195.0, 66.23, 7.4), rgb255(192.0, 62.81, 9.8), rgb255(189.0, 59.38, 12.3), rgb255(186.0, 55.96, 14.8), rgb255(183.0, 52.54, 17.2), rgb255(180.0, 49.12, 19.7), rgb255(177.0, 45.69, 22.2), rgb255(174.0, 42.27, 24.6), rgb255(171.0, 38.85, 27.1), rgb255(168.0, 35.42, 29.5), rgb255(165.0, 32.00, 32.0), rgb255(161.1, 29.54, 32.0), rgb255(157.2, 27.08, 32.0), rgb255(153.2, 24.62, 32.0), rgb255(149.3, 22.15, 32.0), rgb255(145.4, 19.69, 32.0), rgb255(141.5, 17.23, 32.0), rgb255(137.5, 14.77, 32.0), rgb255(133.6, 12.31, 32.0), rgb255(129.7, 9.85, 32.0), rgb255(125.8, 7.38, 32.0), rgb255(121.8, 4.92, 32.0), rgb255(117.9, 2.46, 32.0), rgb255(114.0, 0.00, 32.0), rgb255(115.0, 9.81, 41.3), rgb255(116.1, 19.62, 50.6), rgb255(117.1, 29.42, 59.9), rgb255(118.2, 39.23, 69.2), rgb255(119.2, 49.04, 78.5), rgb255(120.2, 58.85, 87.8), rgb255(121.3, 68.65, 97.2), rgb255(122.3, 78.46, 106.5), rgb255(123.3, 88.27, 115.8), rgb255(124.4, 98.08, 125.1), rgb255(125.4, 107.89, 134.4), rgb255(126.5, 117.69, 143.7), rgb255(127.5, 127.50, 153.0), rgb255(131.4, 131.42, 156.9), rgb255(135.3, 135.35, 160.8), rgb255(139.3, 139.27, 164.8), rgb255(143.2, 143.19, 168.7), rgb255(147.1, 147.12, 172.6), rgb255(151.0, 151.04, 176.5), rgb255(155.0, 154.96, 180.5), rgb255(158.9, 158.88, 184.4), rgb255(162.8, 162.81, 188.3), rgb255(166.7, 166.73, 192.2), rgb255(170.7, 170.65, 196.2), rgb255(174.6, 174.58, 200.1), rgb255(178.5, 178.50, 204.0), rgb255(180.5, 180.46, 206.0), rgb255(182.4, 182.42, 207.9), rgb255(184.4, 184.38, 209.9), rgb255(186.3, 186.35, 211.8), rgb255(188.3, 188.31, 213.8), rgb255(190.3, 190.27, 215.8), rgb255(192.2, 192.23, 217.7), rgb255(194.2, 194.19, 219.7), rgb255(196.2, 196.15, 221.7), rgb255(198.1, 198.12, 223.6), rgb255(200.1, 200.08, 225.6), rgb255(202.0, 202.04, 227.5), rgb255(204.0, 204.00, 229.5), rgb255(206.0, 205.96, 230.5), rgb255(207.9, 207.92, 231.5), rgb255(209.9, 209.88, 232.4), rgb255(211.8, 211.85, 233.4), rgb255(213.8, 213.81, 234.4), rgb255(215.8, 215.77, 235.4), rgb255(217.7, 217.73, 236.4), rgb255(219.7, 219.69, 237.3), rgb255(221.7, 221.65, 238.3), rgb255(223.6, 223.62, 239.3), rgb255(225.6, 225.58, 240.3), rgb255(227.5, 227.54, 241.3), rgb255(229.5, 229.50, 242.2), rgb255(230.5, 230.48, 242.8), rgb255(231.5, 231.46, 243.4), rgb255(232.4, 232.44, 244.0), rgb255(233.4, 233.42, 244.6), rgb255(234.4, 234.40, 245.2), rgb255(235.4, 235.38, 245.8), rgb255(236.4, 236.37, 246.4), rgb255(237.3, 237.35, 247.0), rgb255(238.3, 238.33, 247.5), rgb255(239.3, 239.31, 248.1), rgb255(240.3, 240.29, 248.7), rgb255(241.3, 241.27, 249.3), rgb255(242.2, 242.25, 249.9), rgb255(242.6, 242.64, 250.1), rgb255(243.0, 243.04, 250.3), rgb255(243.4, 243.43, 250.5), rgb255(243.8, 243.82, 250.7), rgb255(244.2, 244.21, 250.9), rgb255(244.6, 244.60, 251.1), rgb255(245.0, 245.00, 251.3), rgb255(245.4, 245.39, 251.5), rgb255(245.8, 245.78, 251.7), rgb255(246.2, 246.17, 251.9), rgb255(246.6, 246.57, 252.1), rgb255(247.0, 246.96, 252.3), rgb255(247.4, 247.35, 252.4), rgb255(247.8, 247.81, 252.7), rgb255(248.3, 248.28, 252.9), rgb255(248.7, 248.74, 253.1), rgb255(249.2, 249.21, 253.4), rgb255(249.7, 249.67, 253.6), rgb255(250.1, 250.13, 253.8), rgb255(250.6, 250.60, 254.1), rgb255(251.1, 251.06, 254.3), rgb255(251.5, 251.52, 254.5), rgb255(252.0, 251.99, 254.8), rgb255(252.4, 252.45, 255.0));
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
  if ( irr != 0 ) {  // irregular points
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
      density = image(x, y, z, Range(zmin, zmax), p);
   else
      density = image(x, y, z, Automatic, p);
  
   if(xlog)
      palette(zlabel, density, (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/20.), ymin), (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/8.), ymax), Right, p, PaletteTicks("")); 
   else
      palette(zlabel, density, (xmax+(xmax-xmin)/20., ymin), (xmax+(xmax-xmin)/8., ymax), Right, p, PaletteTicks("")); 
   return ndata;
  }
  else{ // regular points
   int nxy[] = fin.dimension(2);  // nx, ny
     int nx = nxy[0];
     int ny = nxy[1];
   real[][] z ;
   z = new real[nx][ny];
   z = fin.dimension(nx, ny);
   bounds density;
  if( zmin < zmax )
      density = image( z, Range(zmin, zmax), (xmin,ymin), (xmax, ymax), p);
   else
      density = image( z, Automatic, (xmin,ymin), (xmax, ymax), p);

   if(xlog)
      palette(zlabel, density, (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/20.), ymin), (pow10(log10(xmax)+(log10(xmax)-log10(xmin))/8.), ymax), Right, p, PaletteTicks("")); 
   else
      palette(zlabel, density, (xmax+(xmax-xmin)/20., ymin), (xmax+(xmax-xmin)/8., ymax), Right, p, PaletteTicks(""));
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
    else if(block == "CURVE"){
       nlines = plot_curve(fin);
       write(stdout, 'a curve is plotted from ' + ((string) nlines ) + ' points.\n');}
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
       nlines = plot_expand(fin);}
    else if(block == "LEGEND"){
       nlines = plot_legend(fin);
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
     clip( (xmincoor, ymincoor) -- (xmaxcoor, ymincoor) -- (xmaxcoor, ymaxcoor) -- (xmincoor, ymaxcoor) -- cycle );
  if(caption !=  '')
   label( caption, ( xmincoor*0.5+xmaxcoor*0.5, ymaxcoor+(ymaxcoor-ymincoor)*0.06 ) );
 if(topaxis == 0){
    xaxis(xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax, p=coorpen, ticks=LeftTicks, above=true);
    xaxis("", axis=YEqualsCenter(cymax, false), xmin = cxmin, xmax = cxmax, p=coorpen, ticks=RightTicksNoLabel, above=true);}
else{
    xaxis(xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax,  p=coorpen, ticks=LeftTicks, above=true);
    plot_topaxis();}

if(rightaxis == 0){
   yaxis(ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, p = coorpen, ticks=RightTicks, above=true);
   yaxis("",  axis=XEqualsCenter(cxmax, false), ymin = cymin, ymax = cymax, p=coorpen, ticks=LeftTicksNoLabel, above=true);}
else{ 
   yaxis(ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, p=coorpen, ticks=RightTicks, above=true);
   plot_rightaxis();}

}


void set_scales(){
if(zlog){
 if(xlog && ylog)
   scale(Log, Log, Log);
 else{
    if(xlog)
       scale(Log, Linear, Log);
    else if(ylog)
       scale(Linear, Log, Log);
    else
       scale(Linear, Linear, Log);}}
else{
  if(xlog && ylog)
    scale(Log, Log, Linear);
  else{
    if(xlog)
       scale(Log, Linear, Linear);
    else if(ylog)
       scale(Linear, Log, Linear);
    else
       scale(Linear, Linear, Linear); }}
}

// =============================== Main Routine===============================
//load the file
file fconf =input(name = "asyplot.config", check = false);
textfile = fconf;
if(textfile == "") textfile = getstring(prompt="Enter the 2d image text file: ");
file fin=input(textfile);
// =============================================================================
//read in width and height of the figure, in unit inch
real t[] = fin.dimension(2); 
size(t[0]*inch, t[1]*inch, IgnoreAspect);
// =============================================================================
//read in the captioin, x label, y label
caption = fetch_string(fin);
xlabel = fetch_string(fin);
ylabel = fetch_string(fin);
// =========================================================
//setting xlog, ylog, zlog
int dologs[] = fin.dimension(3); // xlog, ylog, zlog
xlog = (dologs[0] != 0);
ylog = (dologs[1] != 0);
zlog = (dologs[2] != 0);
set_scales();
//===================================
// do clipping?
int i = fin;
doclip = (i != 0);
//==================================================
// read in x, y limits
real[] t;
t = new real [2];
t = fin.dimension(2); //xmin, xmax
cxmin = t[0];
cxmax = t[1];
xmin_adjust = (cxmin >= infty);
xmax_adjust = (cxmax <= -infty);
t = fin.dimension(2); // ymin, ymax
cymin = t[0];
cymax = t[1];
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
int nblocks = fin;
if(nblocks > 0){  //plot the first n blocks
 for(int iblock = 0; iblock < nblocks; ++iblock){
  if(! plot_block(fin)) break;}}
else{ //plot all
  while(plot_block(fin));
 }
// plot the axes
plot_axes();










