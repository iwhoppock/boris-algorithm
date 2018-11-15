# Boris Algorithm or Boris Pusher or Boris Rotation

The Boris algorithm for numerically tracing non-relativistic charged particles in electromagnetic fields. The basic algorithm may be found in either the python or Matlab scripts. The C file allows for multiple particles in a Maxwellian velocity distribution. I have included a matlab plotting script for the C output (.csv format).

## The C Code

To compile the C code, use autotools: (1) autoreconf -i (2) ./configure (3) make (4) ./boris

NB: This is parallelised using OpenMP (shared memory parallel processing). The user may disable the header <omp.h> and the line "#pragma omp parallel for".

## Fortran90

A fortran90 version is available upon request. 
