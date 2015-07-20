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

all: $(LIB) UDG GetPeaks Stack MSMAP ByProd GetInfo GetTheo SST

SWS:	$(OBJ) sample_widths.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) sample_widths.o $(LINK)

MM:	$(OBJ) map_minkowski.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_minkowski.o $(LINK)

LLPW:	$(OBJ) lowl_power.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) lowl_power.o $(LINK)

LLSCAN:	$(OBJ) lowl_scan.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) lowl_scan.o $(LINK)


LLCORR:	$(OBJ) lowl_corr.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) lowl_corr.o $(LINK)

GetTheo:  $(OBJ) gettheo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) gettheo.o $(LINK)

ACTTheo:  $(OBJ) acttheo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) acttheo.o $(LINK)

GetPeaks:  $(OBJ) exportpeaks.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) exportpeaks.o $(LINK)

FindExt: $(OBJ) find_extrema.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) find_extrema.o $(LINK)

GetInfo: $(OBJ) getinfo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getinfo.o $(LINK)

UDG: $(OBJ) udgrade.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) udgrade.o $(LINK)

GetCl: $(OBJ) getcls.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getcls.o $(LINK)

Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

CCards: $(OBJ) change_cards.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) change_cards.o $(LINK)


MPat: $(OBJ) merge_patches.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) merge_patches.o $(LINK)

Inpaint: $(OBJ) testinp.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testinp.o $(LINK)


F2H: $(OBJ) flat2heal.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat2heal.o $(LINK)

Rotate: $(OBJ) rotate_map.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) rotate_map.o $(LINK)

Corr: $(OBJ) corr_map.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) corr_map.o $(LINK)

Blend: $(OBJ) blend_map.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) blend_map.o $(LINK)



SimuMaps: $(OBJ) simulate_maps.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simulate_maps.o $(LINK)

CatP: $(OBJ) peakcat.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) peakcat.o $(LINK)


FullStack: $(OBJ) stack_fullsky.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) stack_fullsky.o $(LINK)


Stack: $(OBJ) dostack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) dostack.o $(LINK)

StackFil: $(OBJ) stack_filaments.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) stack_filaments.o $(LINK)

MST: 	$(OBJ) msst.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) msst.o $(LINK)

SST: $(OBJ) sstack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) sstack.o $(LINK)


ByProd: $(OBJ) map_by_products.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_by_products.o $(LINK)


MSMAP: $(OBJ) map_merge_split.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_merge_split.o $(LINK)

FSpots:  $(OBJ) flat_exportspots.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_exportspots.o $(LINK)

FStack:	 $(OBJ) flat_stack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_stack.o $(LINK)

FSim:	 $(OBJ) flat_sim.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_sim.o $(LINK)

FMerge:	 $(OBJ) flat_merge.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_merge.o $(LINK)

DoShell: $(OBJ) shells.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) shells.o $(LINK)	

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test FMerge Fsim FStack FSpots UDG Stack MSMAP SimuMaps DoShell SST MST ByProd GetInfo GetTheo GetPeaks Rotate Corr Inpaint LLPW



gettheo.o: $(LIB) $(OBJ)

sstack.o: $(LIB) $(OBJ)

test.o: $(LIB) $(OBJ)

testinp.o: $(LIB) $(OBJ)

exportpeaks.o: $(LIB) $(OBJ)

dostack.o : $(LIB) $(OBJ)

flat_exportspots.o: $(LIB) $(OBJ)

flat_stack.o: $(LIB) $(OBJ)

flat_merge.o: $(LIB) $(OBJ)

trimcmp.o: $(LIB) $(OBJ)

map_merge_split.o: $(LIB) $(OBJ)

map_by_products.o: $(LIB) $(OBJ)

corr_map.o:	$(LIB) $(OBJ)

rotate_map.o:	$(LIB) $(OBJ)

healgeo.o: stacking_obj.o

