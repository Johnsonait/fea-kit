#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>

#include "linearsystem.h"
#include "body.h"
#include "quadrature.h"
#include "matrix.h"

void Log(std::string message)
{
	std::cout<< "\n"<< message << "\n";
}

//Base class to define all subsequent elements
class Element
{
private:

public:

};

class TetrahedralElement : public Element
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
			b.push_back(shape_derivatives[chosen_node][i]);

			for (int j = 0; j < nodes[0].size(); j++) //Columns of jacobian matrix
			{
				double sum = 0;
				for (int m = 0; m < nodes.size(); m++) //Number of nodes (x = x1*N1 + ... +xm*Nm)
				{
					sum += nodes[m][i]*shape_derivatives[m][j];
				}
				mat[i].push_back(sum);
			}
		}
		LinearSystem temp_system(mat, b);
		//Update globe_shape_derivatives with solution jacobian equation
		//Can create B matrix after getting global shape derivatives
		temp_system.Solve(global_shape_derivatives[chosen_node]); //Update global shape derivates
	}

	//Constructs 6x12 elemental B matrix
	//Requires global shape function derivatives
	Matrix ConstructBMatrix()
	{
		std::vector<std::vector<double>> Mat;

		Mat.resize(6); //Set rows of Mat to 6
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

	double Jacobian_Det() //Calculate the Jacobian determinant
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

	void GetGlobalShapeDerivatives(double zeta, double eta, double mu)
	{
		for (int m = 0; m < nodes.size(); m++)
		{
			Jacobian(0,0,0,m); //zeta, eta, and mu need to be set to specific values for more complex elements
		}
	}

public:
	TetrahedralElement() //Default constructor
	{

	}
	TetrahedralElement(std::vector<std::vector<double>>& body_nodes, int element[4])
	{
		for (int i = 0; i < body_nodes.size(); i++)
		{
			nodes.push_back(body_nodes[element[i]]);
			global_shape_derivatives.push_back({});
		}

		jacobian_det = Jacobian_Det();
		//Need to call Jacobian to find global shape derivatives before running ConstructBMatrix()
		GetGlobalShapeDerivatives(0,0,0);
		b_matrix = ConstructBMatrix();
		b_matrix.PrintMatrix();
		b_matrix.Transpose();
		b_matrix.PrintMatrix();
	}
};

//Class used to read and parse model data, parameters, and boundary conditions
class Reader 
{
private:
	std::string mesh_file_name;
	std::string solver_input_name;
	std::vector<std::string> instructions;

	int node_numbers;

	void ReadSolverInput()
	{

	}

	void ReadMesh()
	{
		std::ifstream inFile(mesh_file_name);
		std::string in = "";

		int instruction_index = 0;
		std::string instruction;
		while (inFile >> in)
		{
			if (in[0] == '#')
			{
				instructions.push_back(in);
				instruction = in;
			}
			else if (instruction == "#NODES")
			{
				double x, y, z;
				x = std::stod(in);
				inFile >> y;
				inFile >> z;
				body.AddNode(x, y, z);
			}
			else if (instruction == "#ELEMENTS")
			{

			}

			instruction_index++;
		}
		inFile.close();
	}

	void ReadBoundaries()
	{

	}

public:
	Body body;

	Reader(std::string mesh_name,std::string input_name)
	{
		mesh_file_name = mesh_name;
		solver_input_name = input_name;
	}

	void Write()
	{

	}

};

//Class to solve problems in linear elasticity 
class ElasticSolids
{
private:
	std::vector<Element> elements;

public:
	ElasticSolids()//Default constructors
	{

	}
};

//Quick function to test Quadrature class
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

	std::vector<std::vector<double>> nodey = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
	int elem[4] = {0,1,2,3};

	TetrahedralElement elley(nodey, elem);

}

