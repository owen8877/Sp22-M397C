--- Feb 8, 2011, P.G. Martinsson

This directory contains three self-contained files that demonstrate how
the FMM works. Observe that these are stripped down versions that have
been optimized for clarity and not efficiency. 

NOTE: For actual computations, do NOT use these codes.

The files provided are:

--- main_simplefmm
    This is an absolutely minimal code. 
    It computes all operators on the fly, and can only handle the log-kernel.

--- main_simplefmm_visualize
    This is the same code as "main_simplefmm" but it visualizes each step.

--- main_fmm
    This is a code that is a little closer to a "real" FMM.
    It precomputes all the translation operators that are recurring, for instance.
    Still, it is very much an educational code that is not optimized for speed.

--- main_fmm_timing
    This is the same code as "main_fmm" but with some additional lines
    that time the execution of different parts of the code.

======== DATA STRUCTURES ==========

The tree of boxes is encoded in a "cell-object" called NODES.
Each column stores the information associated with one node in the tree.
Specifically, for each box "ibox" it provides the following information:
   NODES{01,ibox} = geometry information for ibox [length,x1_center,x2_center]
   NODES{02,ibox} = the level of ibox
   NODES{03,ibox} = the parent of ibox
   NODES{04,ibox} = list of the children of ibox
   NODES{05,ibox} = number of children
   NODES{06,ibox} = the index list for ibox is ind = NODES{6,ibox}-1+(1:NODES{7,ibox})
   NODES{07,ibox} = the index list for ibox is ind = NODES{6,ibox}-1+(1:NODES{7,ibox})
   NODES{10,ibox} = the list of neighbors of ibox ("list 1")
   NODES{12,ibox} = the interaction list of ibox  ("list 2")
   NODES{14,ibox} = "list 3"
   NODES{16,ibox} = "list 4"
   NODES{40,ibox} = the diagonal block A(J,J) where J = NODES{6,ibox}-1+(1:NODES{7,ibox})
   NODES{41,ibox}{j} = blocks A(J_ibox,J_jbox) for jbox in the neighbor list.
   NODES{43,ibox}{j} = blocks A(J_ibox,J_jbox) for jbox in "list 3".
   NODES{44,ibox}{j} = blocks A(J_ibox,J_jbox) for jbox in "list 4".
   NODES{46,ibox} = T_ofs for ibox
   NODES{47,ibox} = T_tfi for ibox

The geometry of a box is specified by a vector "box_geom" as follows:
   box_geom(1) = sidelength of the box
   box_geom(2) = x1-coordinate of center of box
   box_geom(3) = x2-coordinate of center of box

The data that is _always_ precomputed is stored in the "cell-object" T_OPS.
For each level in the tree, three types of objects are stored:
   T_OPS{1,ilevel}{igeom} = T_ofo   for respective position "igeom"
   T_OPS{2,ilevel}{igeom} = T_ifo   for respective position "igeom"
   T_OPS{3,ilevel}{igeom} = T_ifi   for respective position "igeom"
For T_ofo and T_ifi, igeom is defined as follows
   igeom=2   igeom=4
   igeom=1   igeom=3
For T_ifo, igeom is defined as follows:
    7   14   21   28   35   42   49
    6   13   20   27   34   41   48
    5   12                  40   47
    4   11       ibox       39   46
    3   10                  38   45
    2    9   16   23   30   37   44
    1    8   15   22   29   36   43
In other words, igeom = (i1+3)*7 + i2 + 4 where (i1,i2) are the coordinates
of the box location.

