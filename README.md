COMPOS
======

A tool to compare the trabecular structure of bones across multiple 
specimens

Method
------

We describe the method in 2D for simplicity

1. Given a set of segments describing the contour of bones, triangulate 
   each bone. We expect each triangulation to be unique to each bone due
   to differences in shape.

2. Identify a landmark point on the external contour. Compute a landmark 
   point in the interior, for instance by computing the center of mass of
   the external segments. The external landmark point should correspond to
   the same anatomical region -- this may require human intervention.
   
3. Solve the Laplace equation for two fields: rho and theta. Rho has zero 
   boundary condition on the internal landmark point, and rho = 1 on the 
   external contour. To compute theta, we draw line form the internal 
   landmark to the external landmark point and apply theta = 0 on one 
   side of the line and theta = 1 on the other. Together, rho and theta
   uniquely identify a position on the bone, irrespective of differences
   in sizes, rotation, and elongation. 
   
4. With 1, 2, and 3 in place we can find any location on a bone specimen,
   which matches a location on another, reference bone. From the location 
   on the reference bone compute rho and theta (requires interpolation). 
   Then one finds the x, y coordinates that map to rho and theta on the 
   specimen bone. This requires an iterative method, for instance Newton.
   

   
