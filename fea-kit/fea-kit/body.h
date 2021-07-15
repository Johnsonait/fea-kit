#pragma once

struct Body
{
public:
	std::string name;

	//Array storing the coordinates of nodes, nodes[i] is a vector of 3 double coordinates X Y Z
	//the global nodal number is uniquely determined by the index
	std::vector<std::vector<double>> nodes;

	//Array storing the nodes associated with a number of traction boundaries. Each seperate boundary
	//is tied to the index of the node-storing vector
	std::vector<std::vector<int>> traction_boundary_nodes; //Each entry is a list of nodes, the index ties
	//the list of nodes with the associated traction stored in traction_boundaries
	std::vector<double> traction_boundaries; //Define each traction as [Tx Ty Tz] in array of tractions

	//Array storing the nodes associated with a number of displacement boundaries. Each seperate boundary
	//is tied to the index of the node-storing vector
	std::vector<std::vector<int>> displacement_boundary_nodes; //Same principle as traction_boundary_nodes
	std::vector<double> displacement_boundaries; //Define each displacement [u v w]


	std::vector<std::vector<double>> displacement; //Store the displacement result for each node, []
	std::vector<std::vector<double>> stress; //Store stress values (Voight) [Sxx Syy Szz Szy Szx Sxy]
	std::vector<std::vector<double>> elastic_matrix; //Store the elasticity matrix for the body C
	std::vector<double> temperature;//Store temperatures at nodes

	double elastic_modulus;
	double poisson;

	Body();

	void GetStrain();

	void GetStress();

	void OutputData();

	void AddNode(double x, double y, double z);
};