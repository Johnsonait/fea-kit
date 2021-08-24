#include "triangular_shell.h"



TriangularShell::TriangularShell(std::vector<std::vector<double>>& body_nodes, std::vector<uint32_t>& element)
{
	for (int i = 0; i < body_nodes.size(); i++)
	{
		nodes.push_back(body_nodes[element[i]]);
		global_shape_derivatives.push_back({});
	}

	jacobian_det = JacobianDet();
}

TriangularShell::~TriangularShell() {}

double TriangularShell::JacobianDet()
{
	double x[3] = {0,0,0};
	double y[3] = {0,0,0};
	double z[3] = {0,0,0};

	for (int i = 0; i < nodes.size(); ++i)
	{
		x[i] = nodes[i][0];
		y[i] = nodes[i][1];
		z[i] = nodes[i][2];
	}
	return  ((x[1] - x[0]) * (y[2] - y[0])) - ((y[1] - y[0]) * (x[2] - x[0]));
}

double TriangularShell::ShapeFunction(const double& zeta, const double& eta, const  double& mu, const uint32_t& m)
{
	switch (m)
	{
	case 0:
		return 1 - zeta - eta;
	case 1:
		return zeta;
	case 2:
		return eta;
	default:
		break;
	}
}

double TriangularShell::ShapeFunctionDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m, const uint32_t& direction)
{
	//m refers to one of the 4 possible nodes, 1,2,3,4
	//Direction refers to the 3 local coordinates of zeta, eta,mu
	switch (m)
	{
	//N_0
	case 0:
		switch (direction)
		{
		case 1:
			return -1;
		case 2:
			return -1;
		default:
			break;
		}
	//N_1
	case 1:
		switch (direction)
		{
		case 1:
			return 1;
		case 2:
			return 0;
		default:
			break;
		}
	//N_2
	case 2:
		switch (direction)
		{
		case 1:
			return 0;
		case 2:
			return 1;
		default:
			break;
		}
	default:
		break;
	}
}

void TriangularShell::CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const double& mu, const size_t& node_num)
{

}

const double& TriangularShell::GetJacobianDet(double, double, double)
{
	return JacobianDet();
}