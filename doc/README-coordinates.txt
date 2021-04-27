Last updated: 2016 May 01
------------

Internally, SCID-TDSE works in the local coordinate system, while the observables are
typically needed in the laboratory system. The transformation can be at times confusing,
introducing hard-to-debug errors. Here are the rules:

1. Local coordinate system is chosen such that:
   1a. The local Z axis is along the direction determined by the polar coordinates
       theta and phi.
   1b. The local Y axis remains within the laboratory XY plane

2. This choice of the local coordinate system corresponds to the following sequence of
   coordinate-system rotations:
   2a. Rotate laboratory system around laboratory Z axis by (phi)
   2b. Rotate the new coordinate system around the new Y axis by (theta)

3a. The transformation matrices for transforming a vector given in the laboratory
   coordinate system to the local coordinate system is returned by:

    MathRotationMatrix((/phi,theta,0./),rm)

3b. To transform a vector from laboratory to local coordinate system:

    vec_local = matmul(rm,vec_lab)

3c. To transform a vector from local to laboratory coordinate system:

    vec_lab = matmul(transpose(rm),vec_local)

4. The explicit form of the transformation matrix returned by MathRotationMatrix
   is:

    / Cos theta  0  -Sin theta \ /  Cos phi  Sin phi   0 \
    |    0       1      0      | | -Sin phi  Cos phi   0 |
    \ Sin theta  0   Cos theta / \     0        0      1 /
         
