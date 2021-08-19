#pragma once
#include <iostream>

#include "element.h"

class TetrahedralElement :
	public Element
{
private:
	std::vector<std::vector<double>> jacobian; //Store the 3x3 jacobian matrix
	double jacobian_det; //Store element jacobian value (can potentially change in space)
	//Function that uses jacobian to find global derivatives of shape function 
	//for constructing B matrix
	double JacobianDet() override;

public:
	//Constructors
	TetrahedralElement();
	TetrahedralElement(std::vector<std::vector<double>>& body_nodes,std::vector<uint32_t>&);
	//Destructor
	~TetrahedralElement();

	//Element shape functions in local coordinate system
	//Simple for tets but added for potentential future uses
	double ShapeFunction(const double&, const double&, const  double&, const uint32_t&) override;
	double ShapeFunctionDerivatives(const double&, const double&, const  double&, const uint32_t&, const uint32_t&) override;
	void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const double& mu, const size_t& node_num) override;

	const double& GetJacobianDet(double, double, double) override;

};