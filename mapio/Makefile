SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopmapio.a
OBJ = sphere.o readfits.o fitswrap.o healgeo.o zeta3d.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

all: UDG GetSpots MSMAP Stack ByProd GetInfo


GetInfo: $(OBJ) getinfo.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getinfo.o $(LINK)

UDG: $(OBJ) udgrade.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) udgrade.o $(LINK)


GetCl: $(OBJ) getcls.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getcls.o $(LINK)

TrimC: $(OBJ) trimcmp.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) trimcmp.o $(LINK)

QEB: $(OBJ) quickeb.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) quickeb.o $(LINK)


Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

SimuMaps: $(OBJ) simulate_maps.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simulate_maps.o $(LINK)



GetSpots: $(OBJ) exportspots.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) exportspots.o $(LINK)

Stack: $(OBJ) dostack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) dostack.o $(LINK)


MultStack: $(OBJ) multstack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) multstack.o $(LINK)


HA: $(OBJ) hast.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) hast.o $(LINK)


SST: $(OBJ) sstack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) sstack.o $(LINK)

POSTHA: $(OBJ) post_ha.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) post_ha.o $(LINK)

POSTSST: $(OBJ) post_sst.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) post_sst.o $(LINK)

ByProd: $(OBJ) map_by_products.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_by_products.o $(LINK)


MSMAP: $(OBJ) map_merge_split.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) map_merge_split.o $(LINK)

InpEB: $(OBJ) simeb.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simeb.o $(LINK)

InpT: $(OBJ) simt.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) simt.o $(LINK)


FSpots:  $(OBJ) flat_exportspots.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_exportspots.o $(LINK)

FStack:	 $(OBJ) flat_stack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_stack.o $(LINK)

FSim:	 $(OBJ) flat_sim.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_sim.o $(LINK)

FMerge:	 $(OBJ) flat_merge.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) flat_merge.o $(LINK)

Shells: $(OBJ) shells.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) shells.o $(LINK)	

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test FMerge Fsim FStack FSpots InpEB UDG Stack GetSpots MSMAP InpEB InpT TrimC GetCl QEB MultStack Shells POSTHA HA SST POSTSST ByProd GetInfo


test.o: $(LIB)

exportspots.o: $(LIB)

dostack.o :	$(LIB)

flat_exportspots.o: $(LIB)

flat_stack.o: $(LIB)

asymmetry.o: $(LIB)

ensemble_asymmetry.o: $(LIB)

flat_merge.o: $(LIB)

trimcmp.o: $(LIB)

multstack.o: $(LIB)
