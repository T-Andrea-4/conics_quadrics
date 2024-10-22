A C static library that, given a non-degenerate conic or quadric, finds the canonical form of the section/surface, the center or vertex (if there is one) and the isometry that maps the surface in canonical form. To do that the library implements two functions: conicCanonical and quadricCanonical. 

conicCanonical takes as arguments the matrix of doubles (3 X 3) that represents the conic section (called A), a second matrix of doubles (3 X 3) that after the function call will store the canonical form coefficients (called ACanonical), a third matrix of doubles (3 X 3) that will store the isometry (roto-translation) that maps the conic in canonical form (called isometry), an array of doubles (called centreOrVertex) that will store the coordinates, x and y, of the center (if the conic is an ellipse or a hyperbola) or of the vertex (if the conic is a parabola) of the conic section. The function returns ‘p’, ‘h’, ‘e’, ‘d’ if section is a parabola, a hyperbola, an ellipse, a degenerated conic, respectively. 

quadricCanonical takes as arguments three 4 X 4 matrixes,  (the one that represents the quadric, the one that will store the coefficients of the canonical form and the one that will store the isometry), an array that will store the coordinates (x, y, z) of the center (if the quadric has it), and a pointer to an int variable that, after the function call,  will be1  if the conic has a center, otherwise, 0.

There is another set of functions (implemented in the gauss.c file) that deal with linear algebra and that are used to accomplish the task of the two functions described above.

 Diagonalization.c file contains the implementation of the function used to find the eigenvalues of matrixes using the Jacobi eigenvalues algorithm.
