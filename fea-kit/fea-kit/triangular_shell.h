#pragma once
#include "element.h"

class TriangularShell :
    public Element
{
private:

	double JacobianDet() override;
public:
	TriangularShell() = default;
	TriangularShell(std::vector<std::vector<double>>& body_nodes, std::vector<uint32_t>&);

	~TriangularShell();

	double ShapeFunction(const double&, const double&, const  double&, const uint32_t&) override;
	double ShapeFunctionDerivatives(const double&, const double&, const  double&, const uint32_t&, const uint32_t&) override;
	void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m) override;

	const double& GetJacobianDet(double, double, double) override;
};

