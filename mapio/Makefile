SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
include ../configure.in
include ../compile_rules.in

LIB = libcoopmapio.a
OBJ = sphere.o readfits.o fitswrap.o healgeo.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

all: UDG HSM GetSpots Stack MSMAP InpEB InpT Ena Asym FStack FSim FSpots FMerge

UDG: $(OBJ) udgrade.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) udgrade.o $(LINK)


HSM: $(OBJ) half_sphere_mask.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) half_sphere_mask.o $(LINK)



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

Asym: $(OBJ) asymmetry.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) asymmetry.o $(LINK)

Ena: $(OBJ) ensemble_asymmetry.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) ensemble_asymmetry.o $(LINK)


HSL: $(OBJ) hemisphere_stack_linearfit.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) hemisphere_stack_linearfit.o $(LINK)

POSTHSL: $(OBJ) post_hsl.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) post_hsl.o $(LINK)


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

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test FMerge Fsim FStack FSpots InpEB UDG Stack GetSpots MSMAP HSM InpEB InpT Ena Asym TrimC GetCl QEB MultStack HSL


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
