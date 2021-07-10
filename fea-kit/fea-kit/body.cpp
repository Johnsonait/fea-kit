#include <vector>

#include "body.h"


	//Array storing the coordinates of nodes, nodes[i] is a vector of 3 double coordinates X Y Z
	//the global nodal number is uniquely determined by the index
	std::vector<std::vector<double>> nodes;

	//Array storing the nodes associated with a number of traction boundaries. Each seperate boundary
	//is tied to the index of the node-storing vector
	std::vector<std::vector<int>> traction_boundary_nodes;

	//Array storing the nodes associated with a number of displacement boundaries. Each seperate boundary
	//is tied to the index of the node-storing vector
	std::vector<std::vector<int>> displacement_boundary_nodes;

	std::vector<std::vector<double>> displacement;
	std::vector<std::vector<double>> stress;
	std::vector<std::vector<double>> elastic_matrix;
	std::vector<double> temperature;

	void Body::GetStrain()
	{

	}

	void Body::GetStress()
	{

	}

	void Body::OutputData()
	{

	}
