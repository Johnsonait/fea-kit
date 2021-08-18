#include "tetrahedral_element.h"


//Constructors
TetrahedralElement::TetrahedralElement()
{
	jacobian = {};
	jacobian_det = 0;
}

TetrahedralElement::TetrahedralElement(std::vector<std::vector<double>>& body_nodes, std::vector<uint32_t>& element)
{
	for (int i = 0; i < body_nodes.size(); i++)
	{
		nodes.push_back(body_nodes[element[i]]);
		global_shape_derivatives.push_back({});
	}

	jacobian_det = JacobianDet();
}
//Destructor
TetrahedralElement::~TetrahedralElement() {}


//Useful
double TetrahedralElement::ShapeFunction(const double& zeta,const double& eta,const double& mu,const uint32_t& m)
{
	switch (m)
	{
	case 1:
		return 1 - zeta - eta - mu;
	case 2:
		return zeta;
	case 3:
		return eta;
	case 4:
		return mu;
	default:
		break;
	}
}

double TetrahedralElement::ShapeFunctionDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m, const uint32_t& direction)
{
	//m refers to one of the 4 possible nodes, 1,2,3,4
	//Direction refers to the 3 local coordinates of zeta, eta,mu
	switch (m)
	{
	case 1:
		switch (direction)
		{
		case 1:
			return -1;
		case 2: 
			return -1;
		case 3:
			return -1;
		default:
			break;
		}
	case 2:
		switch (direction)
		{
		case 1:
			return 1;
		case 2:
			return 0;
		case 3:
			return 0;
		default:
			break;
		}
	case 3:
		switch (direction)
		{
		case 1:
			return 0;
		case 2:
			return 1;
		case 3:
			return 0;
		default:
			break;
		}
	case 4:
		switch (direction)
		{
		case 1:
			return 0;
		case 2:
			return 0;
		case 3:
			return 1;
		default:
			break;
		}
	default:
		break;
	}
}

//Function that uses jacobian to find global derivatives of shape function 
//for constructing B matrix
void TetrahedralElement::CalcGlobalShapeDerivatives(const double& zeta,const double& eta,const double& mu, const size_t& node_num)
{
	/*
	for (size_t m = 0; m < node_num; ++m)
	{
		std::vector<std::vector<double>> mat = { {},{},{} }; //Matrix for jacobian system solution
		std::vector<std::vector<double>> b = { {} };//Vector "solution" for jacobian system 

		// For all coordinates associated with nodes (3)
		// Seems iffy, double check this
		for (int i = 0; i < nodes[0].size(); i++) //Rows of jacobian matrix 
		{
			//Update b vector with chosen node shape derivatives
			b[0].push_back(ShapeFunctionDerivatives(0, 0, 0, m, i));

			for (int j = 0; j < nodes[0].size(); j++) //Columns of jacobian matrix
			{
				double sum = 0;
				for (int m = 0; m < nodes.size(); m++) //Number of nodes (x = x1*N1 + ... +xm*Nm)
				{
					sum += nodes[m][i] * ShapeFunctionDerivatives(0, 0, 0, m, j);
				}
				mat[i].push_back(sum);
			}
		}
		LinearSystem temp_system(mat, b);
		//Update globe_shape_derivatives with solution jacobian equation
		//Can create B matrix after getting global shape derivatives
		temp_system.Solve(global_shape_derivatives[m]); //Update global shape derivates
	}
	*/
	//Here the global shape derivatives are found using a hand-calculated solution, bypasses the expensive linear system procedure above
	//Variables x,y,z defined here to help track indices
	int x = 0, y = 1, z = 2;
	//dN_0/dx
	global_shape_derivatives[0][x] = ((nodes[2][1]*nodes[1][2]) -(nodes[3][1]*nodes[1][2]) -(nodes[1][1]*nodes[2][2]) +(nodes[3][1]*nodes[2][2]) +(nodes[1][1]*nodes[3][2]) -(nodes[2][1]*nodes[3][2])) / jacobian_det;
	//dN_0/dy
	global_shape_derivatives[0][y] = ((-1*(nodes[2][x]*nodes[1][z]))+(nodes[3][x]*nodes[1][z])+(nodes[1][x]*nodes[2][z])-(nodes[3][x]*nodes[2][z])-(nodes[1][x]*nodes[3][z])+(nodes[2][x]*nodes[3][z]))/ jacobian_det;
	//dN_0/dz
	global_shape_derivatives[0][z] = ((nodes[2][x]*nodes[1][y])-(nodes[3][x]*nodes[1][y])-(nodes[1][x]*nodes[2][y])+(nodes[3][x]*nodes[2][y])+(nodes[1][x]*nodes[3][y])-(nodes[2][x]*nodes[3][y]))/ jacobian_det;

	//dN_1/dx
	global_shape_derivatives[1][x] = ((-1*(nodes[2][y]*nodes[0][z]))+(nodes[3][y]*nodes[0][z])+(nodes[0][y]*nodes[2][z])-(nodes[3][y]*nodes[2][z])-(nodes[0][y]*nodes[3][z])+(nodes[2][y]*nodes[3][z]))/ jacobian_det;
	//dN_1/dy
	global_shape_derivatives[1][y] = ((nodes[2][x]*nodes[0][z])-(nodes[3][x]*nodes[0][z])-(nodes[0][x]*nodes[2][z])+(nodes[3][x]*nodes[2][z])+(nodes[0][x]*nodes[3][z])-(nodes[2][x]*nodes[3][z])) / jacobian_det;
	//dN_1/dz
	global_shape_derivatives[1][z] = ((-1*(nodes[2][x]*nodes[0][y]))+(nodes[3][x]*nodes[0][y])+(nodes[0][x]*nodes[2][y])-(nodes[3][x]*nodes[2][y])-(nodes[0][x]*nodes[3][y])+(nodes[2][x]*nodes[3][y])) / jacobian_det;
	
	//dN_2/dx
	global_shape_derivatives[2][x] = ((nodes[1][y]*nodes[0][z])-(nodes[3][y]*nodes[0][z])-(nodes[0][y]*nodes[1][z])+(nodes[3][y]*nodes[1][z])+(nodes[0][y]*nodes[3][z])-(nodes[1][y]*nodes[3][z])) / jacobian_det;
	//dN_2/dy
	global_shape_derivatives[2][y] = ((-1*(nodes[1][x]*nodes[0][z]))+(nodes[3][x]*nodes[0][z])+(nodes[0][x] * nodes[1][z])-(nodes[3][x] * nodes[1][z])-(nodes[0][x] * nodes[3][z])+(nodes[1][x] * nodes[3][x]))/ jacobian_det;
	//dN_2/dz
	global_shape_derivatives[2][z] = ((nodes[1][x]*nodes[0][y])-(nodes[3][x] * nodes[0][y])-(nodes[0][x] * nodes[1][y])+(nodes[3][x] * nodes[1][y])+(nodes[0][x] * nodes[3][y])-(nodes[1][x] * nodes[3][y]))/ jacobian_det;

	//dN_3/dx
	global_shape_derivatives[3][x] = ((-1*(nodes[1][y]*nodes[0][z]))+(nodes[2][y] * nodes[0][z])+(nodes[0][y] * nodes[1][z])-(nodes[2][y] * nodes[1][z])-(nodes[0][y] * nodes[2][z])+(nodes[1][y] * nodes[2][z]))/ jacobian_det;
	//dN_3/dy
	global_shape_derivatives[3][y] = ((nodes[1][x]*nodes[0][z])-(nodes[2][x]*nodes[0][z])-(nodes[0][x] * nodes[1][z])+(nodes[2][x] * nodes[1][z])+(nodes[0][x] * nodes[2][z])-(nodes[1][x] * nodes[2][z])) / jacobian_det;
	//dN_3/dz
	global_shape_derivatives[3][z] = ((-1*(nodes[1][x]*nodes[0][y]))+(nodes[2][x] * nodes[0][y])+(nodes[0][x] * nodes[1][y])-(nodes[2][x] * nodes[1][y])-(nodes[0][x] * nodes[2][y])+(nodes[1][x] * nodes[2][y])) / jacobian_det;
}

double TetrahedralElement::JacobianDet() //Calculate the Jacobian determinant
{
	//Not necessary but makes the determinant calculation much easier to read
	double x[4] = {0,0,0,0};
	double y[4] = {0,0,0,0};
	double z[4] = {0,0,0,0};

	for (int i = 0; i < nodes.size(); ++i)
	{
		x[i] = nodes[i][0];
		y[i] = nodes[i][1];
		z[i] = nodes[i][2];
	}

	return  ((x[1] - x[0]) * (((y[2] - y[0]) * (z[3] - z[0])) - ((y[3] - y[0]) * (z[2] - z[0]))))
		- ((y[1] - y[0]) * (((x[2] - x[0]) * (z[3] - z[0])) - ((x[3] - x[0]) * (z[1] - z[0]))))
		+ ((z[1] - z[0]) * (((x[2] - x[0]) * (y[3] - y[0])) - ((x[3] - x[0]) * (y[2] - y[0]))));
}
