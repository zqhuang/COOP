SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopmapio.a
OBJ = sphere.o stacking_obj.o readfits.o healgeo.o fitswrap.o  zeta3d.o

default: all

$(LIB): $(OBJ)
	ar -r $@ $?

all: $(LIB) UDG GetPeaks Stack MSMAP ByProd GetInfo GetTheo TE2Z

GAUSS: $(OBJ) gaussianize.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) gaussianize.o $(LINK)

PLOTLINES: $(OBJ) plot_lines.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) plot_lines.o $(LINK)

MM0:	$(OBJ) minkowski0.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) minkowski0.o $(LINK)

MM1:	$(OBJ) minkowski1.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) minkowski1.o $(LINK)

GetTheo:  $(OBJ) gettheo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) gettheo.o $(LINK)

ACTTheo:  $(OBJ) acttheo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) acttheo.o $(LINK)

GetPeaks:  $(OBJ) exportpeaks.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) exportpeaks.o $(LINK)

GetInfo: $(OBJ) getinfo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getinfo.o $(LINK)

UDG: $(OBJ) udgrade.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) udgrade.o $(LINK)

Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

CutPatch: $(OBJ) healpix2patch.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) healpix2patch.o $(LINK)

SimuMaps: $(OBJ) simulate_maps.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simulate_maps.o $(LINK)

Stack: $(OBJ) dostack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) dostack.o $(LINK)

FStack: $(OBJ) flat_stack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_stack.o $(LINK)


ByProd: $(OBJ) map_by_products.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_by_products.o $(LINK)

MSMAP: $(OBJ) map_merge_split.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_merge_split.o $(LINK)

TE2Z: $(OBJ) te2zeta.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) te2zeta.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test PLOTLINES MSMAP ByProd Stack SimuMaps Geteaks GetInfo ACTTheo GetTheo MM0 MM1 TE2Z



gettheo.o: $(LIB) $(OBJ)


test.o: $(LIB) $(OBJ)


exportpeaks.o: $(LIB) $(OBJ)

dostack.o : $(LIB) $(OBJ)

map_merge_split.o: $(LIB) $(OBJ)

map_by_products.o: $(LIB) $(OBJ)

rotate_map.o:	$(LIB) $(OBJ)

healgeo.o: stacking_obj.o

te2z.o: $(LIB) $(OBJ)
