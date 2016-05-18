Below is a descriptions of the results files.

Tri-Seq-p#.txt: uniform refinement using provided (python script)
triangular meshes, and full-order (triangular) approximation space

Quad-Seq-p#.txt: uniform refinement using provided (python script)
quadrilateral meshes, and tensor-product (quad) approximation space

QuadTri-Seq-p#.txt: uniform refinement using provided (python script)
quadrilateral meshes, and full-order (triangular) approximation space

HangNode-p*.txt: hanging-node fixed fraction (15%) quadrilateral
refinement starting with ref0 of the provided quad mesh sequence.
Adaptation on drag (bottom wall) using output adjoints.

AdaptP.txt: order (p) refinement, output-based on drag, fixed fraction
(20%), starting from ref1 of the provided quadrilateral mesh sequence.

Tri-MOESS-p#.txt: output-based triangular mesh optimization using
MOESS-like approach with several iterations at roughly fixed dof, and
five dof levels.








