# Boris Algorithm or Boris Pusher or Boris Rotation

The Boris algorithm for numerically tracing non-relativistic charged particles in electromagnetic fields. The basic algorithm may be found in either the python or Matlab scripts. The C file allows for multiple particles in a Maxwellian velocity distribution. I have included a matlab plotting script for the C output (.csv format). I have included a function file in the Matlab format. 

## The C Code

To compile the C code, use autotools: (1) autoreconf -i (2) ./configure (3) make (4) ./boris

NB: This is parallelised using OpenMP (shared memory parallel processing). The user may disable the header <omp.h> and the line "#pragma omp parallel for".

## Fortran90

A fortran90 version is available upon request. 

## Additional Resources

For additional resources for particle pushing, I recommend these two delightful blog posts and their associated references. I believe there is also sample code written in Java: 

- https://www.particleincell.com/2011/velocity-integration/
- https://www.particleincell.com/2011/vxb-rotation/

Dear Ian,

I'm tui from Hiroshima Univeristy, Japan. Now I'm trying to simulate particle under magnetic and electric field by using Boris method for particle pusher. I also applied your code into fortran language and constrcuted the main program. My fortran program assumed particle moving in 1D spatial space with 3D velocity. The magnetic field and electric field is the same direction in z-axis. Anyway, the velocity profiles from my simulation provided flucmmtuated to kinetic energy when I squared those velcoities. Could you give any comments to my program? The main program is available from https://drive.google.com/drive/folders/1vxzcUSz337SQlb5np-PJOYH8C9PLIG58?usp=drive_link . Thank you very much   
