#pragma once

#include "element.h"

class TetrahedralElement :
	public Element
{
private:
	std::vector<std::vector<double>> jacobian; //Storing value of 3x3 jacobian matrix mapping between local and global coordinates
	double jacobian_det;
	//Function that uses jacobian to find global derivatives of shape function 
	//for constructing B matrix
	void Jacobian(double zeta, double eta, double mu, int m);
	double JacobianDet() override;


public:
	//Constructors
	TetrahedralElement();
	TetrahedralElement(std::vector<std::vector<double>>& body_nodes, int element[4]);

	//Element shape functions in local coordinate system
	//Simple for tets but added for potentential future uses
	double ShapeFunction(double, double, double, uint32_t);

	//Accessors
	const std::vector<std::vector<double>>& GetGlobalShapeDerivatives(double zeta, double eta, double mu);
	const std::vector<std::vector<double>>& GetNodes();

	//Mutators
	void AddNode(const std::vector<double>& n);

};