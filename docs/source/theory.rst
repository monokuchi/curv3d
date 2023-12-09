======
Theory
======
This documents the mathematical theory and computational optimziations behind the
rasterization process. Rasterization is the process of drawing pixels on the screen
given some 3D surface(s). In our case, we use Bézier surfaces to substitute the
conventional triangle.

In rasterization, we can loosely split the process up into
a few steps:

1. Camera Projection
2. Pixel Intersection
3. Depth Handling


Camera Projection
=================

To deal with the camera,
