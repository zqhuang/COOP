SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopmapio.a
OBJ = stacking_obj.o healutils.o fitswrap.o  zeta3d.o

default: all

$(LIB): $(OBJ)
	ar -r $@ $?

all: $(LIB) UDG GetPeaks Stack MSMAP GetInfo GetTheo XFASTER ClFormat

flat: FCoadd FMerge FSmooth FSmoothQU FStat FGetPeaks FStack FQU2EB FEB2QU FMask FSim F2Heal Heal2F FPTOP MCoadd FTrim

ITM:  $(OBJ) calc_intermittency.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) calc_intermittency.o $(LINK)


XFASTER: $(OBJ) xfaster.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) xfaster.o $(LINK)

MAP2CLS: $(OBJ) getcls.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getcls.o $(LINK)

GAUSS: $(OBJ) gaussianize.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) gaussianize.o $(LINK)

NAREAS: $(OBJ) count_areas.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) count_areas.o $(LINK)

FDUST: $(OBJ) filterdust.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) filterdust.o $(LINK)

MDUST: $(OBJ) multfildust.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) multfildust.o $(LINK)

INP:	$(OBJ) inpaint.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) inpaint.o $(LINK)

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

ClFormat: $(OBJ) clform.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) clform.o $(LINK)

SimuMaps: $(OBJ) simulate_maps.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simulate_maps.o $(LINK)

Stack: $(OBJ) dostack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) dostack.o $(LINK)

FStack: $(OBJ) flat_stack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_stack.o $(LINK)

FGetPeaks: $(OBJ) flat_exportspots.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_exportspots.o $(LINK)

FStat:  $(OBJ) flat_stat.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_stat.o $(LINK)

FQU2EB:  $(OBJ) flat_qu2eb.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_qu2eb.o $(LINK)

FEB2QU:  $(OBJ) flat_eb2qu.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_eb2qu.o $(LINK)

F2Heal:  $(OBJ) flat2heal.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat2heal.o $(LINK)


Heal2F:  $(OBJ) heal2flat.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) heal2flat.o $(LINK)

FSmooth:  $(OBJ) flat_smooth.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_smooth.o $(LINK)

FSmoothQU:  $(OBJ) flat_smooth_qu.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_smooth_qu.o $(LINK)

FTrim:  $(OBJ) flat_cut_radec.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_cut_radec.o $(LINK)

FMerge:  $(OBJ) flat_merge.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_merge.o $(LINK)

FCoadd:  $(OBJ) flat_coadd.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_coadd.o $(LINK)

F1D:  $(OBJ) flat_1ddistr.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_1ddistr.o $(LINK)

FPTOP:  $(OBJ) flat_ptop.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_ptop.o $(LINK)

FMask:  $(OBJ) flat_mask.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_mask.o $(LINK)

FGenMask:  $(OBJ) flat_genmask.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_genmask.o $(LINK)


FSim:  $(OBJ) flat_sim.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_sim.o $(LINK)

MCoadd:  $(OBJ) matcoadd.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) matcoadd.o $(LINK)


ByProd: $(OBJ) map_by_products.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_by_products.o $(LINK)

MSMAP: $(OBJ) map_merge_split.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_merge_split.o $(LINK)

TE2Z: $(OBJ) te2zeta.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) te2zeta.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test UDG GetPeaks Stack MSMAP ByProd GetInfo GetTheo TE2Z GAUSS MM0 MM1 FStack FStat F2Heal FCoadd ACTTheo FMask FGetPeaks FQU2EB FEB2QU FSmooth FSmoothQU FMerge FGenMask Heal2F FPTOP MCoadd F1D FDUST INP MAP2CLS GAUSS MDUST ClFormat XFASTER ITM SimuMaps



gettheo.o: $(LIB) $(OBJ)


test.o: $(LIB) $(OBJ)


exportpeaks.o: $(LIB) $(OBJ)

dostack.o : $(LIB) $(OBJ)

map_merge_split.o: $(LIB) $(OBJ)

map_by_products.o: $(LIB) $(OBJ)

rotate_map.o:	$(LIB) $(OBJ)

healutils.o: stacking_obj.o

te2z.o: $(LIB) $(OBJ)

