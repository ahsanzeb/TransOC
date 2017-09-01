
objects = basisstates.o main.o hamiltonian.o lists.o modmain.o memory.o \
 hoppings.o dhops.o kinds.o multiply.o annihilation.o kappa.o gamma.o \
 maps.o init.o amplitudes.o Aorder.o 

include make.inc

### execname= $(HOME)/xphotonon
execname= transoc

switch= -O3 -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
# -g -fbounds-check
# -g 
# -fbounds-check
# -static

#-----------------------------------------------------

transoc: $(objects)
	$(FC) -o $(execname) $(switch) $(objects) $(LIB)

#-----------------------------------------------------
# modules
#-----------------------------------------------------
modmain.mod:	modmain.o modmain.f
	$(FC) -c $(switch) modmain.f

lists.mod: lists.o lists.f
	$(FC) -c $(switch) lists.f

kinds.mod: kinds.o kinds.f
	$(FC) -c $(switch) kinds.f

basisstates.mod: modmain.mod basisstates.o basisstates.f
	$(FC) -c $(switch) basisstates.f

hamiltonian.mod: modmain.mod lists.mod basisstates.mod hamiltonian.o hamiltonian.f
	$(FC) -c $(switch)  hamiltonian.f 

dhops.mod: modmain.mod dhops.o dhops.f
	$(FC) -c $(switch) dhops.f

Aorder.mod: Aorder.o Aorder.f
	$(FC) -c $(switch) Aorder.f


annihilation.mod: modmain.mod annihilation.o annihilation.f
	$(FC) -c $(switch) annihilation.f

creation.mod: modmain.mod Aorder.mod creation.o creation.f
	$(FC) -c $(switch) creation.f

kappa.mod: modmain.mod kappa.o kappa.f
	$(FC) -c $(switch) kappa.f

gamma.mod: modmain.mod gamma.o gamma.f
	$(FC) -c $(switch) gamma.f

hoppings.mod: modmain.mod dhops.mod annihilation.mod creation.mod\
 kappa.mod gamma.mod hoppings.o hoppings.f
	$(FC) -c $(switch) hoppings.f

#-----------------------------------------------------
# objects
#-----------------------------------------------------

modmain.o:	modmain.f
	$(FC) -c $(switch) modmain.f


Aorder.o:	Aorder.f
	$(FC) -c $(switch) Aorder.f

lists.o: lists.f
	$(FC) -c $(switch) lists.f

memory.o: memory.f
	$(FC) -c $(switch) memory.f

kinds.o: kinds.f
	$(FC) -c $(switch) kinds.f

multiply.o: multiply.f
	$(FC) -c $(switch) multiply.f


init.o: modmain.mod init.f
	$(FC) -c $(switch) init.f


basisstates.o: modmain.mod basisstates.f
	$(FC) -c $(switch) basisstates.f

maps.o:	modmain.mod maps.f
	$(FC) -c $(switch) maps.f


dhops.o: modmain.mod basisstates.mod lists.mod amplitudes.o dhops.f
	$(FC) -c $(switch) dhops.f

annihilation.o: modmain.mod basisstates.mod amplitudes.o annihilation.f
	$(FC) -c $(switch) annihilation.f

creation.o: modmain.mod basisstates.mod amplitudes.o Aorder.o creation.f
	$(FC) -c $(switch) creation.f

kappa.o: modmain.mod basisstates.mod amplitudes.o kappa.f
	$(FC) -c $(switch) kappa.f

gamma.o: modmain.mod basisstates.mod lists.mod amplitudes.o gamma.f
	$(FC) -c $(switch) gamma.f

main.o:	modmain.mod basisstates.mod hamiltonian.mod lists.mod hoppings.mod main.f
	$(FC) -c $(switch) main.f 

hamiltonian.o: modmain.mod basisstates.mod hamiltonian.f
	$(FC) -c $(switch) hamiltonian.f

amplitudes.o: modmain.mod amplitudes.f
	$(FC) -c $(switch) amplitudes.f

hoppings.o: modmain.mod dhops.mod annihilation.mod creation.mod\
 kappa.mod gamma.mod hoppings.f
	$(FC) -c $(switch) hoppings.f

#-----------------------------------------------------


%.o: %.f

clean:
	rm $(execname) *.o  *.mod
# =================================================================================
# End of the makefile
# =================================================================================

