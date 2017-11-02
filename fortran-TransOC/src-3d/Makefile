
objects = basisstates.o main.o hamiltonian.o lists.o modmain.o memory.o \
 hoppings.o dhops.o modmultiply.o annihilation.o losses.o \
 maps.o init.o amplitudes.o creation.o selection.o \
 rates.o modways.o readinput.o diag.o timestamp.o contacts.o lindblad.o



include make.inc


switch= -O3 -g -fbounds-check -fbacktrace 
#switch= -O3 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal
# -O -Wall -fcheck=all -g -fbacktrace
# -g -fbounds-check
# -g 
# -fbounds-check
# -static

#-----------------------------------------------------

transoc: $(objects)
	$(FC) -o transoc $(switch) $(objects) $(LIB)
	cp transoc ~/bin/.

LIB = libarpack_mac.a -llapack 

#-----------------------------------------------------
# modules
#-----------------------------------------------------
modmain.mod:	modmain.o modmain.f
	$(FC) -c $(switch) modmain.f

readinput.mod: modmain.mod readinput.o readinput.f
	$(FC) -c $(switch) readinput.f

maps.mod: modmain.mod lists.mod maps.o maps.f
	$(FC) -c $(switch) maps.f

diag.mod: modmain.mod diag.f
	$(FC) -c $(switch) diag.f

lists.mod: modmain.mod lists.o lists.f
	$(FC) -c $(switch) lists.f

basisstates.mod: modmain.mod lists.mod basisstates.o basisstates.f
	$(FC) -c $(switch) basisstates.f

hamiltonian.mod: modmain.mod lists.mod basisstates.mod hamiltonian.o hamiltonian.f
	$(FC) -c $(switch)  hamiltonian.f 


lindblad.mod: modmain.mod lindblad.o lindblad.f
	$(FC) -c $(switch)  lindblad.f 

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

amplitudes.mod: modmain.mod basisstates.mod modmultiply.mod amplitudes.f
	$(FC) -c $(switch) amplitudes.f

modmultiply.mod: modmain.mod modmultiply.o modmultiply.f
	$(FC) -c $(switch) modmultiply.f

contacts.mod: modmain.mod basisstates.mod lists.mod amplitudes.mod lists.mod contacts.f
	$(FC) -c $(switch) contacts.f

#-----------------------------------------------------
# objects
#-----------------------------------------------------

modmain.o:	modmain.f
	$(FC) -c $(switch) modmain.f

lists.o: modmain.mod lists.f
	$(FC) -c $(switch) lists.f

readinput.o: modmain.mod readinput.f
	$(FC) -c $(switch) readinput.f

memory.o: memory.f
	$(FC) -c $(switch) memory.f

modmultiply.o: modmultiply.f
	$(FC) -c $(switch) modmultiply.f

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

creation.o: modmain.mod basisstates.mod lists.mod amplitudes.o creation.f
	$(FC) -c $(switch) creation.f


contacts.o: modmain.mod basisstates.mod lists.mod amplitudes.o lists.o contacts.f
	$(FC) -c $(switch) contacts.f

losses.o: modmain.mod basisstates.mod lists.mod amplitudes.o losses.f
	$(FC) -c $(switch) losses.f

lindblad.o: modmain.mod lindblad.f
	$(FC) -c $(switch) lindblad.f

main.o:	modmain.mod basisstates.mod hamiltonian.mod lists.mod \
hoppings.mod init.mod selection.mod rates.mod modways.mod\
 readinput.mod maps.mod diag.mod lindblad.mod main.f 
	$(FC) -c $(switch) main.f 

hamiltonian.o: modmain.mod basisstates.mod hamiltonian.f
	$(FC) -c $(switch) hamiltonian.f

amplitudes.o: modmain.mod modmultiply.mod amplitudes.f
	$(FC) -c $(switch) amplitudes.f

hoppings.o: modmain.mod dhops.mod annihilation.mod creation.mod\
losses.mod contacts.mod hoppings.f
	$(FC) -c $(switch) hoppings.f

rates.o: modmain.mod rates.f
	$(FC) -c $(switch) rates.f

modways.o: modmain.mod lists.mod modways.f
	$(FC) -c $(switch) modways.f

diag.o: modmain.mod diag.f
	$(FC) -c $(switch) diag.f

timestamp.o: timestamp.f
	$(FC) -c $(switch) timestamp.f

#-----------------------------------------------------


%.o: %.f

clean:
	rm transoc *.o  *.mod

# =================================================================================
# End of the makefile
# =================================================================================

