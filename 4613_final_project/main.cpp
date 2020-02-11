
#include <cstring>
#include <iostream>
#include <Eigen/Dense>

#include "readmsh.h"
#include "shapefunctions.h"
#include "inhomogeneousf.h"

int main(int argc, char* argv[])
{
	// the name of the mesh file, copy it in case we need to manipulate it
	char *fname = new char[strlen(argv[1]) + 1];
	strcpy(fname, argv[1]);

	// the name of the boundary file, in case we need to manipulate it
	char *bname = new char[strlen(argv[2]) + 1];
	strcpy(bname, argv[2]);

	/* 
	 * read in the mesh
	 */
	msh m(fname, bname);

	/*
	 * set up the inhomogeneous function, it is not based on user input!
	 */
	inh::constant(0);
	inh::quadratic(0, 0, 3, 0);

	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(m.node_len, m.node_len);
	Eigen::VectorXd f = Eigen::VectorXd::Zero(m.node_len);

	/*
	 * based on the elements that have been read in, we calculate the force vector
	 * and the element stiffness matrix
	 */
	for (int i = 0; i < m.element_len; ++i)
	{
		auto &element = m.elements[i];
		if (element.type == TRIANGLE_ID)
		{
			shape_function{ m.nodes }.triangle(
				element.node_ids[0], element.node_ids[1], element.node_ids[2], 
				K.data(), m.node_len);
		}
	}
	for (int i = 0; i < m.node_len; ++i)
	{
		f(i) = inh::constant();
	}

	Eigen::MatrixXd K_unknowns = Eigen::MatrixXd::Zero(m.node_types.inode_len, m.node_types.inode_len);
	Eigen::VectorXd f_unknowns = Eigen::VectorXd::Zero(m.node_types.inode_len);

	// build each row in smaller stiffness matrix for the unknown quantities
	// based on what we know on the boundaries (or other known points)
	for (int i = 0; i < m.node_types.inode_len; ++i)
	{
		size_t index_i = m.node_types.interior_nodes[i];

		// iterate over all the boundary nodes to sum the known quantities and bring
		// the values over onto the right side with the f vector
		double sum = f(index_i);
		for (int j = 0; j < m.node_types.bnode_len; ++j)
		{
			size_t index_j = m.node_types.boundary_nodes[j];
			double value = m.node_types.boundary_values[j];

			sum -= K(index_i, index_j) * value;
		}
		f_unknowns(i) = sum;

		// iterate over all the unknown quantities to copy the coefficients in the columns
		for (int j = 0; j < m.node_types.inode_len; ++j)
		{
			size_t index_j = m.node_types.interior_nodes[j];
			K_unknowns(i, j) = K(index_i, index_j);
		}
	}


	Eigen::VectorXd u_unknowns = K_unknowns.inverse() * f_unknowns;

	// now that we have the unknown values, build the full vector of heat values
	double *u = new double[m.node_len];
	for (int i = 0; i < m.node_types.bnode_len; ++i)
	{
		size_t index = m.node_types.boundary_nodes[i];
		u[index] = m.node_types.boundary_values[i];
	}
	for (int i = 0; i < m.node_types.inode_len; ++i)
	{
		size_t index = m.node_types.interior_nodes[i];
		u[index] = u_unknowns(i);
	}

	FILE *out;
	if ((out = fopen("results.txt", "w")) == 0)
	{
		printf("results file could not be opened with name '%s'\n", fname);

	}

	fprintf(out, "$map << EOD\n");
	for (int i = 0; i < m.element_len; ++i)
	{
		auto &element = m.elements[i];
		if (element.type == TRIANGLE_ID)
		{
			for (int ii = 0; ii <= 3; ++ii)
			{
				size_t index_i = element.node_ids[ii % 3];
				fprintf(out, "%.2lf %.2lf %.2lf\n", m.nodes[index_i].x, m.nodes[index_i].y, u[index_i]);
			}
		}
		fprintf(out, "\n\n");
	}
	char texname[100];
	int i = 0;
	for (char *c = fname; *c != '.'; ++c)
	{
		texname[i++] = *c;
	}
	texname[i] = '\0';

	fprintf(out, "EOD\n");
	fprintf(out, R"~(

set term epslatex size 5.5,4
set output "%s.tex"

unset key
set xlabel "$x$"
set ylabel "$y$"
set zlabel "heat ($u$)" rotate by 90
set title "Result of Heat Problem"

splot $map with lines

unset output
)~", texname);
}

