#pragma once

#include "readmsh.h"

/*
 * this struct has knowledge of the mesh and uses it to apply the shape functions.
 * the shape functions were previously derived and are specific to the element types
 * calling the functions returns a corresponding element in the LOCAL stiffness matrix
 */
struct shape_function
{
	msh::point *nodes;

	/*
	 * the shape function to get the stiffness matrix entry when the element is a triangle
	 * this function will add the values from the local stiffness matrix to the global stiffness matrix
	 * the first three parameters, 'i', 'j', 'k' are the node indicies from the node list that
	 * correspond to the nodes of the triangular element for which the stiffness matrix is being computed
	 * the parameter 'K' is the 1d array containing stiffness element matrix entries, aligned by row
	 * the 'row_len' parameter is the length of one row in the stiffness matrix (also the number of nodes)
	 */
	inline void triangle(size_t i, size_t j, size_t k, double *K, size_t row_len)
	{
		// record the x and y positions of the nodes in the order they are
		// passed, for convenient access
		double
			x[] = { nodes[i].x, nodes[j].x, nodes[k].x },
			y[] = { nodes[i].y, nodes[j].y, nodes[k].y };

		double a[] = {
			x[1] - x[2],
			x[2] - x[0],
			x[0] - x[1] };

		double b[] = {
			y[1] - y[2],
			y[2] - y[0],
			y[0] - y[1] };

		// compute the jacobian, which is always positive
		double J = abs(a[1] * b[2] - a[2] * b[1]);

		size_t ind[] = { i, j, k };
		for (int ii = 0; ii < 3; ++ii)
		{
			for (int jj = 0; jj < 3; ++jj)
			{
				K[ind[ii] * row_len + ind[jj]] += 1.0 / (2 * J) * (a[ii] * a[jj] + b[ii] * b[jj]);
			}
		}
	}
};


