#include "tetrahedral_element.h"





double TetrahedralElement::ShapeFunction(double zeta, double eta, double mu,uint32_t n)
{
	switch (n)
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

//Function that uses jacobian to find global derivatives of shape function 
//for constructing B matrix
void TetrahedralElement::Jacobian(double zeta, double eta, double mu, int chosen_node)
{
	std::vector<std::vector<double>> mat = { {},{},{} }; //Matrix for jacobian system solution
	std::vector<double> b = {};//Vector "solution" for jacobian system 

	// For all coordinates associated with nodes (3)
	// Seems iffy, double check this
	for (int i = 0; i < nodes[0].size(); i++) //Rows of jacobian matrix 
	{
		//Update b vector with chosen node shape derivatives
		b.push_back(shape_derivatives[chosen_node][i]);

		for (int j = 0; j < nodes[0].size(); j++) //Columns of jacobian matrix
		{
			double sum = 0;
			for (int m = 0; m < nodes.size(); m++) //Number of nodes (x = x1*N1 + ... +xm*Nm)
			{
				sum += nodes[m][i] * shape_derivatives[m][j];
			}
			mat[i].push_back(sum);
		}
	}
	LinearSystem temp_system(mat, b);
	//Update globe_shape_derivatives with solution jacobian equation
	//Can create B matrix after getting global shape derivatives
	temp_system.Solve(global_shape_derivatives[chosen_node]); //Update global shape derivates
}

double TetrahedralElement::Jacobian_Det() //Calculate the Jacobian determinant
{
	//Not necessary but makes the determinant calculation much easier to read
	double x[4] = {};
	double y[4] = {};
	double z[4] = {};

	for (int i = 0; i < nodes[0].size(); i++)
	{
		x[i] = nodes[i][0];
		y[i] = nodes[i][1];
		z[i] = nodes[i][2];

	}

	return  ((x[1] - x[0]) * (((y[2] - y[0]) * (z[3] - z[0])) - ((y[3] - y[0]) * (z[2] - z[0]))))
		- ((y[1] - y[0]) * (((x[2] - x[0]) * (z[3] - z[0])) - ((x[3] - x[0]) * (z[1] - z[0]))))
		+ ((z[1] - z[0]) * (((x[2] - x[0]) * (y[3] - y[0])) - ((x[3] - x[0]) * (y[2] - y[0]))));
}


TetrahedralElement::TetrahedralElement() = default;

TetrahedralElement::TetrahedralElement(std::vector<std::vector<double>>& body_nodes, int element[4])
{
	for (int i = 0; i < body_nodes.size(); i++)
	{
		nodes.push_back(body_nodes[element[i]]);
		global_shape_derivatives.push_back({});
	}

	jacobian_det = Jacobian_Det();
	//Need to call Jacobian to find global shape derivatives before running ConstructBMatrix()
	GetGlobalShapeDerivatives(0, 0, 0);
}

const std::vector<std::vector<double>>& TetrahedralElement::GetGlobalShapeDerivatives(double zeta, double eta, double mu)
{
	for (int m = 0; m < nodes.size(); m++)
	{
		Jacobian(0, 0, 0, m); //zeta, eta, and mu need to be set to specific values for more complex elements
	}
	return global_shape_derivatives;
}

const std::vector<std::vector<double>>& TetrahedralElement::GetNodes()
{
	return nodes;
}