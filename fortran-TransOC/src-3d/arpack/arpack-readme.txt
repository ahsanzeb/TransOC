issue the following command
        zcat arpack96.tar.Z | tar -xvf -
        zcat patch.tar.Z | tar -xvf - 
which will create a directory named ARPACK. Change to the
ARPACK directory and read the file README for further instruction.



edit ARmake.inc

fix issue with function etime in seconds.f in util, if there is one.

run:
make lib

to to examples/sym

dsdrv1.f is the double precision standard eigensolver example routine



