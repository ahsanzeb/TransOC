
objects = basisstates.o main.o hamiltonian.o lists.o modmain.o

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

basisstates.mod: modmain.mod basisstates.o basisstates.f
	$(FC) -c $(switch) basisstates.f

hamiltonian.mod: modmain.mod lists.mod basisstates.mod hamiltonian.o hamiltonian.f
	$(FC) -c $(switch)  hamiltonian.f 


#-----------------------------------------------------
# objects
#-----------------------------------------------------

modmain.o:	modmain.f
	$(FC) -c $(switch) modmain.f

lists.o: lists.f
	$(FC) -c $(switch) lists.f

basisstates.o: modmain.mod basisstates.f
	$(FC) -c $(switch) basisstates.f

main.o:	modmain.mod basisstates.mod hamiltonian.mod lists.mod main.f
	$(FC) -c $(switch) main.f 

hamiltonian.o: modmain.mod basisstates.mod hamiltonian.f
	$(FC) -c $(switch) hamiltonian.f

#-----------------------------------------------------


%.o: %.f

clean:
	rm $(execname) *.o  *.mod
# =================================================================================
# End of the makefile
# =================================================================================

