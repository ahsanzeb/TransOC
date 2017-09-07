
objects = basisstates.o main.o hamiltonian.o lists.o modmain.o memory.o \
 hoppings.o dhops.o kinds.o multiply.o annihilation.o losses.o \
 maps.o init.o amplitudes.o creation.o selection.o rates.o modways.o readinput.o diag.o



include make.inc

execname= transoc

#switch= -O3 -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
switch= -O3
# -O -Wall -fcheck=all -g -fbacktrace
# -g -fbounds-check
# -g 
# -fbounds-check
# -static

#-----------------------------------------------------

transoc: $(objects)
	$(FC) -o $(execname) $(switch) $(objects) $(LIB)

LIB = libarpack_mac.a

#-----------------------------------------------------
# modules
#-----------------------------------------------------
modmain.mod:	modmain.o modmain.f
	$(FC) -c $(switch) modmain.f

readinput.mod: modmain.mod readinput.o readinput.f
	$(FC) -c $(switch) readinput.f

maps.mod: modmain.mod lists.mod maps.f
	$(FC) -c $(switch) maps.f

diag.mod: modmain.mod diag.f
	$(FC) -c $(switch) diag.f

lists.mod: lists.o lists.f
	$(FC) -c $(switch) lists.f

kinds.mod: kinds.o kinds.f
	$(FC) -c $(switch) kinds.f

basisstates.mod: modmain.mod lists.mod basisstates.o basisstates.f
	$(FC) -c $(switch) basisstates.f

hamiltonian.mod: modmain.mod lists.mod basisstates.mod hamiltonian.o hamiltonian.f
	$(FC) -c $(switch)  hamiltonian.f 

dhops.mod: modmain.mod dhops.o dhops.f
	$(FC) -c $(switch) dhops.f

init.mod: modmain.mod maps.mod modways.mod init.o init.f
	$(FC) -c $(switch) init.f


annihilation.mod: modmain.mod annihilation.o annihilation.f
	$(FC) -c $(switch) annihilation.f

creation.mod: modmain.mod creation.o creation.f
	$(FC) -c $(switch) creation.f

losses.mod: modmain.mod losses.o losses.f
	$(FC) -c $(switch) losses.f

hoppings.mod: modmain.mod dhops.mod annihilation.mod creation.mod\
 losses.mod hoppings.o hoppings.f
	$(FC) -c $(switch) hoppings.f

rates.mod: modmain.mod rates.o rates.f
	$(FC) -c $(switch) rates.f

selection.mod: modmain.mod selection.f
	$(FC) -c $(switch) selection.f

modways.mod: modmain.mod lists.mod modways.o modways.f
	$(FC) -c $(switch) modways.f

#-----------------------------------------------------
# objects
#-----------------------------------------------------

modmain.o:	modmain.f
	$(FC) -c $(switch) modmain.f

lists.o: lists.f
	$(FC) -c $(switch) lists.f

readinput.o: modmain.mod readinput.f
	$(FC) -c $(switch) readinput.f

memory.o: memory.f
	$(FC) -c $(switch) memory.f

kinds.o: kinds.f
	$(FC) -c $(switch) kinds.f

multiply.o: multiply.f
	$(FC) -c $(switch) multiply.f


init.o: modmain.mod maps.mod modways.mod init.f
	$(FC) -c $(switch) init.f

basisstates.o: modmain.mod lists.mod basisstates.f
	$(FC) -c $(switch) basisstates.f

selection.o: modmain.mod selection.f
	$(FC) -c $(switch) selection.f

maps.o:	modmain.mod lists.mod maps.f
	$(FC) -c $(switch) maps.f

dhops.o: modmain.mod basisstates.mod lists.mod amplitudes.o dhops.f
	$(FC) -c $(switch) dhops.f

annihilation.o: modmain.mod basisstates.mod amplitudes.o annihilation.f
	$(FC) -c $(switch) annihilation.f

creation.o: modmain.mod basisstates.mod amplitudes.o creation.f
	$(FC) -c $(switch) creation.f


losses.o: modmain.mod basisstates.mod lists.mod amplitudes.o losses.f
	$(FC) -c $(switch) losses.f

main.o:	modmain.mod basisstates.mod hamiltonian.mod lists.mod \
hoppings.mod init.mod selection.mod rates.mod modways.mod\
 readinput.mod maps.mod diag.mod main.f 
	$(FC) -c $(switch) main.f 

hamiltonian.o: modmain.mod basisstates.mod hamiltonian.f
	$(FC) -c $(switch) hamiltonian.f

amplitudes.o: modmain.mod amplitudes.f
	$(FC) -c $(switch) amplitudes.f

hoppings.o: modmain.mod dhops.mod annihilation.mod creation.mod\
losses.mod hoppings.f
	$(FC) -c $(switch) hoppings.f

rates.o: modmain.mod rates.f
	$(FC) -c $(switch) rates.f

modways.o: modmain.mod lists.mod modways.f
	$(FC) -c $(switch) modways.f

diag.o: modmain.mod diag.f
	$(FC) -c $(switch) diag.f

#-----------------------------------------------------


%.o: %.f

clean:
	rm $(execname) *.o  *.mod
# =================================================================================
# End of the makefile
# =================================================================================

