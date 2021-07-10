#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>

#include "linearsystem.h"
#include "body.h"

void Log(std::string message)
{
	std::cout<< "\n"<< message << "\n";
}

const std::vector<std::vector<double>> QUADRATURE_POINTS = { 
	{0},
	{-0.5773502691896257,0.5773502691896257},
	{0,-0.7745966692414834,0.7745966692414834},
	{0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538},
};

const std::vector<std::vector<double>> QUADRATURE_WEIGHTS = { 
	{2},
	{1,1},
	{0.8888888888888888,0.5555555555555556,0.5555555555555556},
	{-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526},
};

//Class to implement Gaussian Quadrature on an arbitrary function that takes and returns doubles
class Quadrature
{
public:
	//Gaussian Quadrature for one-dimensional integration
	double IntegrateOneD(const int& points, std::function<double (double)> Func)
	{
		int index = points - 1;
		double result = 0.0;
		for (int i = 0; i < points; i++)
		{
			result += QUADRATURE_WEIGHTS[index][i] * Func(QUADRATURE_POINTS[index][i]);
		}
		return result;
	}
	//Gaussian Quadrature for two-dimensional integration
	double IntegrateTwoD(const int& points, std::function<double (double,double)> Func)
	{
		int index = points - 1;
		double result = 0.0;
		for (int i = 0; i < points; i++)
		{
			for (int j = 0; j< points; j++)
			{
				result += QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j] 
					* Func(QUADRATURE_POINTS[index][i],QUADRATURE_POINTS[index][j]);
			}
		}
		return result;
	}
	//Gaussian Quadrature for three-dimensional integration
	double IntegrateThreeD(const int& points, std::function<double (double,double,double)> Func)
	{
		int index = points - 1;
		double result = 0.0;
		for (int i = 0; i < points; i++)
		{
			for (int j = 0; j < points; j++)
			{
				for (int k = 0; k < points; k++)
				{
					result += QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j] * QUADRATURE_WEIGHTS[index][k]
						* Func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j], QUADRATURE_POINTS[index][k]);
				}
			}
		}
		return result;
	}
};

class Matrix
{
private:
	std::vector<std::vector<double>> matrix;

	int CountRows()
	{
		return matrix.size();
	}

	int CountCols()
	{
		return matrix[0].size();
	}

public:
	Matrix()
	{
	}

	Matrix(std::vector<std::vector<double>>& mat)
	{
		matrix = mat;
	}
	Matrix operator * (Matrix& B)
	{
		std::vector<std::vector<double>> temp;

		for (int i = 0; i < CountRows(); i++)
		{
			std::vector<double> row;
			for (int j = 0; j < B.CountCols(); j++)
			{
				double sum = 0;
				for (int k = 0; k < CountCols(); k++)
				{
					sum += matrix[i][k] * B.matrix[k][j];
				}
				row.push_back(sum);
			}
			temp.push_back(row);
		}
		Matrix ret(temp);
		return ret;
	}

	void Transpose()
	{
		std::vector<std::vector<double>> temp;
		temp.resize(CountCols());

		for (int i = 0; i < CountRows(); i++)
		{
			for (int j = 0; j < CountCols(); j++)
			{
				temp[j].push_back(matrix[i][j]);
			}
		}
		matrix = temp;
	}

	void PrintMatrix()
	{
		std::cout << "\n";
		for (int i = 0; i < matrix.size(); i++)
		{
			for (int j = 0; j < matrix[0].size(); j++)
			{
				std::cout << matrix[i][j] << " ";
			}
			std::cout << "\n";
		}
	}

};


class TetrahedralElement 
{
private:
	std::vector<std::vector<double>> nodes; //Array storing node coordinates as X Y Z
	std::vector<std::vector<double>> jacobian; //Storing value of 3x3 jacobian matrix mapping between local and global coordinates
	double jacobian_det;
	std::vector<std::vector<double>> shape_derivatives = { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} }; //Storing shape derivates (dN1/dzeta etc) 4 nodes, 3 directions gives 12 values X Y Z to zeta eta mu
	std::vector<std::vector<double>> global_shape_derivatives; //Store global shape derivatives for each node
	Matrix b_matrix; //Storing B matrix of element, tetrahedral element has 6x12 B matrix

	//Element shape functions in local coordinate system
	//Simple for tets but added for potentential future uses
	double N1(double zeta, double eta, double mu) { return 1 - zeta - eta - mu; }
	double N2(double zeta) { return zeta; }
	double N3(double eta) { return eta; }
	double N4(double mu) { return mu; }
	
	//Function that uses jacobian to find global derivatives of shape function 
	//for constructing B matrix
	void Jacobian(double zeta, double eta, double mu, int chosen_node) 
	{
		std::vector<std::vector<double>> mat = { {},{},{} }; //Matrix for jacobian system solution
		std::vector<double> b = {};//Vector "solution" for jacobian system 

		// For all coordinates associated with nodes (3)
		// Seems iffy, double check this
		for (int i = 0; i < nodes[0].size(); i++) //Rows of jacobian matrix 
		{
			//Update b vector with chosen node shape derivatives
			b.push_back(shape_derivatives[chosen_node-1][i]);

			for (int j = 0; j < nodes[0].size(); j++) //Columns of jacobian matrix
			{
				double sum = 0;
				for (int m = 0; m < nodes.size(); m++) //Number of nodes (x = x1*N1 + ... +xm*Nm)
				{
					sum += nodes[m][i]*shape_derivatives[m][j];
				}
				mat[i][j] = sum;
			}
		}
		LinearSystem temp_system(mat, b);
		//Update globe_shape_derivatives with solution jacobian equation
		//Can create B matrix after getting global shape derivatives
		temp_system.Solve(global_shape_derivatives[chosen_node-1]);
	}

	//Constructs 6x12 elemental B matrix
	//Requires global shape function derivatives
	Matrix ConstructBMatrix()
	{
		std::vector<std::vector<double>> Mat;
		Mat.resize(6);
		for (int m = 0; m < nodes.size(); m++)
		{
			for (int i = 0; i < 6; i++) //6 rows
			{
				std::vector<std::vector<double>> sub_matrix = {
					{global_shape_derivatives[m][0],0,0},
					{0,global_shape_derivatives[m][1],0},
					{0,0,global_shape_derivatives[m][2]},
					{0,global_shape_derivatives[m][2],global_shape_derivatives[m][1]},
					{global_shape_derivatives[m][2],0,global_shape_derivatives[m][0]},
					{global_shape_derivatives[m][1],global_shape_derivatives[m][0],0}
				};
				for (int j = 0; j < 3; j++) //For col of sub-matrix
				{
					Mat[i].push_back(sub_matrix[i][j]); //Add rows to full B matrix
				}
			}
		}
		Matrix constructed_matrix(Mat);
		return constructed_matrix;
	}

	double Jacobian_Det()
	{
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

public:
	TetrahedralElement(std::vector<std::vector<double>>& body_nodes, int element[4])
	{
		for (int i = 0; i < body_nodes.size(); i++)
		{
			nodes.push_back(body_nodes[element[i]]);
			global_shape_derivatives.push_back({});
		}

		jacobian_det = Jacobian_Det();
		b_matrix = ConstructBMatrix();
	}
};

class Reader 
{

};

double Funky(double x, double y, double z)
{
	return x * x * y * y * z * z;
}

int main()
{
	/*
	std::vector<std::vector<double>> Mat = {{1,0,0,0},{2,4,3,0.5},{1,3,5,2},{1,1,1,1}};
	std::vector<double> b = {1,1,1,1};

	std::vector<double> x;
	 
	LinearSystem System(Mat, b);
	System.Solve(x);
	System.PrintSol(x);
	*/

	//Quadrature integrator;

	//std::cout << integrator.IntegrateThreeD(2,Funky);

	std::vector<std::vector<double>> init = { {1,1,1},{2,2,2},{1,1,1},{1,3,1},{0,1,0} };
	std::vector<std::vector<double>> doubtit = { {1},{1},{1} };
	Matrix A(init);
	Matrix B(doubtit);

	A.PrintMatrix();
	B.PrintMatrix();

	Matrix C = A*B;
	C.PrintMatrix();
	C.Transpose();
	C.PrintMatrix();

}
