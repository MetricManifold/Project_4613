#pragma once

#include <iostream>
#include <cstring>
#include <algorithm>

#define BUFFER 256

#define POINT_ID 15
#define EDGE_ID 1
#define TRIANGLE_ID 2
#define QUADRANGLE_ID 3
#define TRIANGLE_2ND_ID 9
#define QUADRANGLE_2ND_ID 16

struct msh
{
	/*
	 * a type definition of a node, it allows me to easily set the properties
	 * of a node and access those properties (the position in 2d space)
	 */
	using point = struct { double x, y; };

	/* 
	 * a type definition for the element, it allows me to easily set the nodes indices
	 * that the element is composed of, and specify its type
	 */
	using element = struct { size_t node_ids[8]; int type; };

	point *nodes;						// a list of all the nodes in the system
	element *elements;					// a list of all the elements in the system
	size_t node_len, element_len;		// the length of the node and element lists

	/*
	 * this structure is to conveniently separate the full list of nodes into interior nodes
	 * and boundary nodes (based on an input file), then conveniently access the elements
	 * this object name is prefixed with underscores because it shouldn't be used outside this class
	 */
	struct __node_types
	{
		size_t *boundary_nodes;			// the indices of the nodes on the boundary
		size_t *interior_nodes;			// the indices of the interior nodes in the system
		double *boundary_values;		// the values of u that are specified on the boundary
		size_t bnode_len, inode_len;	// the length of the array containing boundary and interior nodes

		/*
		 * this function populates the index mapping for the boundary and interior nodes
		 * the parameter 'fname' is the name of the file containining the boundary value specifications
		 * the parameter 'nodes' is the full list of nodes in the system
		 * 'node_len' is the length of the node list
		 */
		void set_node_types(const char *fname, point *nodes, size_t node_len)
		{
			FILE *f;
			if ((f = fopen(fname, "r")) == 0)
			{
				printf("boundary file does not exist with name '%s'\n", fname);
			}

			char line[BUFFER];

			// the first line contains the number of nodes that are on the boundary
			fgets(line, BUFFER, f);
			sscanf(line, "%zd\n", &bnode_len);
			inode_len = node_len - bnode_len;

			// create the array for all the nodes
			interior_nodes = new size_t[inode_len];
			boundary_nodes = new size_t[bnode_len];
			boundary_values = new double[bnode_len];

			// parse all the boundary points specified in this file
			size_t index = 0;
			while (fgets(line, BUFFER, f))
			{
				size_t node_id;
				double node_value;
				sscanf(line, "%zd %lf\n", &node_id, &node_value);

				// ensure to normalize the indexing to C++ index values
				boundary_nodes[index] = node_id - 1;
				boundary_values[index] = node_value;
				++index;
			}

			// copy the boundary nodes into a list that we will sort. Sorting this list means
			// we can 'fill in' the gaps as the interior nodes
			size_t *boundary_nodes_sort = new size_t[bnode_len];
			std::copy(boundary_nodes, boundary_nodes + bnode_len, boundary_nodes_sort);
			std::sort(boundary_nodes_sort, boundary_nodes_sort + bnode_len);

			index = 0;
			int start = 0, end = (int)boundary_nodes_sort[0];
			for (int i = 1; i < bnode_len + 1; ++i)
			{
				// from one boundary node to the next boundary node, fill in the in the node indices
				// for the interior nodes
				for (int j = start; j < end; ++j)
				{
					interior_nodes[index++] = j;
				}
				start = end + 1;
				end = (i == bnode_len) ? (int)node_len : (int)boundary_nodes_sort[i];
			}
			// fill in everything after the boundary nodes
			for (int j = start; j < end; ++j)
			{
				interior_nodes[index++] = j;
			}

			delete[] boundary_nodes_sort;
		}

		~__node_types() 
		{ 
			delete[] boundary_nodes; 
			delete[] boundary_values;
			delete[] interior_nodes;
		}
	} node_types;


	~msh() { delete[] nodes; delete[] elements; }

	/* 
	 * the constructor to this object. By passing the name of the mesh file, as well as the name
	 * of the file containing the boundary specifications, the information is parsed and the arrays
	 * containing the nodes are populated.
	 * the 'fname' parameter is the name of the mesh file
	 * the 'boundary_fname' parameter is the name of the boundary specification file
	 */
	msh(const char *fname, const char *boundary_fname) : node_len{ 0 }, element_len{ 0 }
	{
		FILE *f;
		if ((f = fopen(fname, "r")) == 0)
		{
			printf("mesh file does not exist with name '%s'\n", fname);
		}

		char line[BUFFER];
		while (fgets(line, BUFFER, f))
		{
			// in the mesh file, the node section has the title "$Nodes", so we will parse
			// all the nodes while in this section
			if (strcmp(line, "$Nodes\n") == 0)
			{
				fgets(line, BUFFER, f);
				sscanf(line, "%*d %zd %*d %*d\n", &node_len);
				nodes = new point[node_len];

				// keep parsing until we get to the end of the nodes section
				while (fgets(line, BUFFER, f), strcmp(line, "$EndNodes\n") != 0)
				{
					int c;

					// read in the entity parameters
					sscanf(line, "%*d %*d %*d %d\n", &c);

					// check first if there are any elements to read undrer this entity
					if (c > 0)
					{

						size_t *indices = new size_t[c];
						size_t *indices_it = indices;

						// the file format first gives the node ids, read these first
						for (int i = 0; i < c; ++i)
						{
							fgets(line, BUFFER, f);
							sscanf(line, "%zd\n", indices_it++);
						}

						// next we read the coordinates, and insert in the same order as node ids
						indices_it = indices;
						for (int i = 0; i < c; ++i)
						{
							double x, y;

							fgets(line, BUFFER, f);
							sscanf(line, "%lf %lf %*lf\n", &x, &y);

							// it is important to insert at ONE LESS than the given size, because 
							// the indexing for gmsh starts at 1, but for C++ it is 0
							nodes[*indices_it++ - 1] = { x, y };
						}

						delete[] indices;
					}
				}
			}

			// this is the elements section which contains our edge bounded shapes
			if (strcmp(line, "$Elements\n") == 0)
			{
				fgets(line, BUFFER, f);

				sscanf(line, "%*d %zd %*d %*d\n", &element_len);
				elements = new element[element_len];
				size_t index = 0;

				// keep reading this section until we get to the end
				while (fgets(line, BUFFER, f), strcmp(line, "$EndElements\n") != 0)
				{
					int type, c;
					sscanf(line, "%*d %*d %d %d\n", &type, &c);

					for (int i = 0; i < c; ++i)
					{
						fgets(line, BUFFER, f);
						if (type == POINT_ID) // single point
						{
							// we don't do anything for a point element
							// this is an explicit point added by a user in gmsh to 
							// create the system, but is already taken from the list of nodes
						}
						else if (type == EDGE_ID) // an edge
						{
							// we don't do anything for an edge element
							// an edge is added when forming the mesh in gmsh, but doesn't
							// contribute to the solution other than bounding an element
						}
						else if (type == TRIANGLE_ID) // a triangle
						{
							sscanf(line, "%*zd %zd %zd %zd\n",
								elements[index].node_ids,
								elements[index].node_ids + 1,
								elements[index].node_ids + 2);
							elements[index].type = type;

							// due to indexing, make sure to increment each point mapping back once
							elements[index].node_ids[0]--;
							elements[index].node_ids[1]--;
							elements[index].node_ids[2]--;

							++index;
						}
						else if (type == QUADRANGLE_ID) // a quadrangle
						{
							sscanf(line, "%*zd %zd %zd %zd %zd\n",
								elements[index].node_ids,
								elements[index].node_ids + 1,
								elements[index].node_ids + 2,
								elements[index].node_ids + 3);
							elements[index].type = type;

							// due to indexing, make sure to increment each point mapping back once
							elements[index].node_ids[0]--;
							elements[index].node_ids[1]--;
							elements[index].node_ids[2]--;
							elements[index].node_ids[3]--;

							++index;
						}
						else if (type == TRIANGLE_2ND_ID) // a 6 node 2nd order triangle
						{
							sscanf(line, "%*zd %zd %zd %zd %zd %zd %zd\n",
								elements[index].node_ids,
								elements[index].node_ids + 1,
								elements[index].node_ids + 2,
								elements[index].node_ids + 3,
								elements[index].node_ids + 4,
								elements[index].node_ids + 5);
							elements[index].type = type;

							// due to indexing, make sure to increment each point mapping back once
							elements[index].node_ids[0]--;
							elements[index].node_ids[1]--;
							elements[index].node_ids[2]--;
							elements[index].node_ids[3]--;
							elements[index].node_ids[4]--;
							elements[index].node_ids[5]--;

							++index;
						}
						else if (type == QUADRANGLE_2ND_ID) // an 8 node 2nd order quadrangle
						{
							sscanf(line, "%*zd %zd %zd %zd %zd %zd %zd %zd %zd\n",
								elements[index].node_ids,
								elements[index].node_ids + 1,
								elements[index].node_ids + 2,
								elements[index].node_ids + 3,
								elements[index].node_ids + 4,
								elements[index].node_ids + 5,
								elements[index].node_ids + 6,
								elements[index].node_ids + 7);
							elements[index].type = type;

							// due to indexing, make sure to increment each point mapping back once
							elements[index].node_ids[0]--;
							elements[index].node_ids[1]--;
							elements[index].node_ids[2]--;
							elements[index].node_ids[3]--;
							elements[index].node_ids[4]--;
							elements[index].node_ids[5]--;
							elements[index].node_ids[6]--;
							elements[index].node_ids[7]--;

							++index;
						}
					}
				}
				element_len = index;
			}
		}

		/*
		 * now that we've parsed all the nodes from the msh file, we can establish the boundary
		 * and interior nodes
		 */
		node_types.set_node_types(boundary_fname, nodes, node_len);
	}
};
