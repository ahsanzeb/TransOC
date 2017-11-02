
perl -pi -w -e 's/write\(std/if(std)write\(\*/g;' *.f


perl -pi -w -e 's/if\(std\)write/write/g;' *.f

perl -pi -w -e 's/use modmain, only: std/	/g;' *.f

perl -pi -w -e 's/std,/ /g;' *.f

perl -pi -w -e 's/,std/ /g;' *.f


