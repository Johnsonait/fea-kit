#pragma once
#include "element.h"

class TetrahedralElement :
	protected Element
{
private:
	std::vector<std::vector<double>> jacobian; //Storing value of 3x3 jacobian matrix mapping between local and global coordinates
	double jacobian_det;
	std::vector<std::vector<double>> shape_derivatives = { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} }; //Storing shape derivates (dN1/dzeta etc) 4 nodes, 3 directions gives 12 values X Y Z to zeta eta mu
	std::vector<std::vector<double>> global_shape_derivatives; //Store global shape derivatives for each node
	std::vector<std::vector<double>> nodes; //Array storing node coordinates as X Y Z

	//Element shape functions in local coordinate system
	//Simple for tets but added for potentential future uses
	double ShapeFunction(double,double,double,uint32_t);

	//Function that uses jacobian to find global derivatives of shape function 
	//for constructing B matrix
	void Jacobian(double zeta, double eta, double mu, int chosen_node);

	double Jacobian_Det();


public:

	TetrahedralElement() = default;

	TetrahedralElement(std::vector<std::vector<double>>& body_nodes, int element[4]);

	const std::vector<std::vector<double>>& GetGlobalShapeDerivatives(double zeta, double eta, double mu);

	const std::vector<std::vector<double>>& GetNodes();
};