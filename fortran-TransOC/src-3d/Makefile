##################################
#					makefile for tranosc3d
##################################

FC 	= mpif90
FFLAGS= -O3 -g -fbounds-check -fbacktrace
LOADLIBES = libarpack_mac.a -llapack 

SRC = basisstates.f main.f hamiltonian.f lists.f modmain.f memory.f \
 hoppings.f dhops.f modmultiply.f annihilation.f losses.f \
 maps.f init.f amplitudes.f creation.f selection.f \
 rates.f modways.f readinput.f diag.f timestamp.f contacts.f coulomb.f

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

transoc3d: $(OBJ)
	$(FC) -o transoc3d $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp transoc3d ~/bin/.


# Here are the compile steps
 
modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 

lists.o: lists.f modmain.o 
	 $(FC) -c $(FFLAGS) lists.f 

maps.o: maps.f lists.o modmain.o
	 $(FC) -c $(FFLAGS) maps.f 

modways.o: modways.f modmain.o lists.o 
	 $(FC) -c $(FFLAGS) modways.f 

init.o: init.f modmain.o maps.o modways.o 
	 $(FC) -c $(FFLAGS) init.f 

basisstates.o: basisstates.f modmain.o lists.o 
	 $(FC) -c $(FFLAGS) basisstates.f 

hamiltonian.o: hamiltonian.f modmain.o lists.o basisstates.o
	 $(FC) -c $(FFLAGS) hamiltonian.f 

modmultiply.o: modmultiply.f  
	 $(FC) -c $(FFLAGS) modmultiply.f 

amplitudes.o: amplitudes.f modmultiply.o modmain.o
	 $(FC) -c $(FFLAGS) amplitudes.f 

dhops.o: dhops.f amplitudes.o lists.o modmain.o basisstates.o 
	 $(FC) -c $(FFLAGS) dhops.f 

creation.o: creation.f amplitudes.o lists.o basisstates.o
	 $(FC) -c $(FFLAGS) creation.f 

annihilation.o: annihilation.f amplitudes.o modmain.o basisstates.o
	 $(FC) -c $(FFLAGS) annihilation.f 

losses.o: losses.f amplitudes.o modmain.o basisstates.o lists.o 
	 $(FC) -c $(FFLAGS) losses.f 

contacts.o: contacts.f amplitudes.o lists.o basisstates.o modmain.o
	 $(FC) -c $(FFLAGS) contacts.f 

hoppings.o: hoppings.f modmain.o dhops.o creation.o annihilation.o losses.o contacts.o 
	 $(FC) -c $(FFLAGS) hoppings.f 

rates.o: rates.f modmain.o 
	 $(FC) -c $(FFLAGS) rates.f 

selection.o: selection.f modmain.o 
	 $(FC) -c $(FFLAGS) selection.f 

readinput.o: readinput.f modmain.o 
	 $(FC) -c $(FFLAGS) readinput.f 

diag.o: diag.f modmain.o
	 $(FC) -c $(FFLAGS) diag.f 

coulomb.o: coulomb.f modmain.o 
	 $(FC) -c $(FFLAGS) coulomb.f 

main.o: main.f modmain.o init.o basisstates.o hamiltonian.o hoppings.o rates.o selection.o modways.o readinput.o maps.o diag.o coulomb.o
	 $(FC) -c $(FFLAGS) main.f 

timestamp.o: timestamp.f
	 $(FC) -c $(FFLAGS) timestamp.f 
	

clean:
	rm -rf *.o *.mod transoc3d
  
