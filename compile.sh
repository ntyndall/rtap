echo "Compiling rtap: run with ./rtap.x"
gfortran routines/rtap.f90 \
         routines/waveorder.f90 \
      -o rtap.x

