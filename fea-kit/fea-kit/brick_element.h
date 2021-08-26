#pragma once

#include <functional>

#include "element.h"

class LinearElasticSolids;

class BrickElement :
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
	BrickElement();
	BrickElement(const std::vector<std::vector<double>>&, std::vector<uint32_t>&);
	//Destructor
	~BrickElement();

	//Element shape functions in local coordinate system
	//Simple for tets but added for potentential future uses
	double ShapeFunction(const double&, const double&, const  double&, const uint32_t&) override;
	double ShapeFunctionDerivatives(const double&, const double&, const  double&, const uint32_t&, const uint32_t&) override;
	void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const double& mu, const size_t& node_num) override;

	const double& GetJacobianDet(double, double, double) override;
	Matrix& Integrate(const int& points, std::function<Matrix& (double, double, double, std::shared_ptr<Element>, LinearElasticSolids*)> func, const Matrix& mat, std::shared_ptr<Element>, LinearElasticSolids*) override;
};

