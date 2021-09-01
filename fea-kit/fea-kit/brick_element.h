#pragma once

#include <functional>

#include "element.h"

class LinearElasticSolids;

class BrickElement :
    public Element
{
private:
	std::vector<std::vector<double>> jacobian; //Store the 3x3 jacobian matrix

	std::vector<std::vector<double>>& Jacobian(const double&, const double&, const double&) override;
	void Jacobian(const double& xsi, const double& eta, const double& zeta, std::vector<std::vector<double>>& ret);

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
	void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const double& mu) override;

	const std::vector<std::vector<double>>& GetJacobian(const double& xsi, const double& eta, const double& zeta) override;

	//const double& GetJacobianDet() override;

	const double& GetJacobianDet(const double& xsi, const double& eta, const double& zeta) override;
};

