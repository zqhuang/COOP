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
Each block can be either DOTS, LINES, CURVE, LABEL, CONTOUR, CLIP, LEGEND, EXTRA_AXIS, or DENSITY
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
coordinates of dot #1 (2 real numbers)
label #1  (string)
coordinates of dot #2 (2 real numbers)
label #2  (string)
....
coordinates of dot #n (2 real numbers)
label #n  (string)
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

// =============================================================================
//plot labels
int plot_legend(file fin){
  string cstr;
  cstr = fetch_string(fin);
  if(trim_string(cstr) !=""){
     int cols = fin;
     if(trim_string(cstr) == "N")
        add(legend(cols), point(N), 20N, UnFill); 
     else if(trim_string(cstr) == "S")
        add(legend(cols), point(S), 20S, UnFill); 
     else if(trim_string(cstr) == "W")
        add(legend(cols), point(W), 20W, UnFill);
     else
        add(legend(cols), point(E), 20E, UnFill);}
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
        yaxis(pic, rightaxis_label, XEqualsRight(cxmax), black, LeftTicks("", begin=false,end=false));
        });}
  else{
    q=secondaryY(new void(picture pic) {
         if(xlog)
            scale(pic,Log,Linear);
         else
            scale(pic,Linear,Linear);
         ylimits(pic, rightaxis_ymin, rightaxis_ymax);
         yaxis(pic, rightaxis_label, XEqualsRight(cxmax), black, LeftTicks("", begin=false,end=false));
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
    xaxis(pic, topaxis_label, YEqualsTop(cymax), black, RightTicks("", begin=false,end=false));
        });}
 else{
   q=secondaryX(new void(picture pic) {
   if(ylog)
       scale(pic, Linear, Log);
   else
      scale(pic, Linear, Linear);
   xlimits(pic, topaxis_xmin, topaxis_xmax);
   xaxis(pic, topaxis_label, YEqualsTop(cymax), black, RightTicks("", begin=false, end=false));
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
  else if(ctbl == "MyColorTable")
    p = Gradient(256, purple, lightblue, mediumblue,  blue, cyan, green, yellow, orange, red, mediumred, lightred);       
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
    xaxis(xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax, ticks=LeftTicks, above=true);
    xaxis("", axis=YEqualsCenter(cymax, false), xmin = cxmin, xmax = cxmax, ticks=RightTicksNoLabel, above=true);}
else{
    xaxis(xlabel, axis=YEqualsCenter(cymin, false), xmin = cxmin, xmax = cxmax, ticks=LeftTicks, above=true);
    plot_topaxis();}

if(rightaxis == 0){
   yaxis(ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, ticks=RightTicks, above=true);
   yaxis("",  axis=XEqualsCenter(cxmax, false), ymin = cymin, ymax = cymax, ticks=LeftTicksNoLabel, above=true);}
else{ 
   yaxis(ylabel,  axis=XEqualsCenter(cxmin, false), ymin = cymin, ymax = cymax, ticks=RightTicks, above=true);
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










