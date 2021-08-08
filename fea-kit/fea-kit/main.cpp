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
#include "reader.h"
#include "element.h"

void Log(std::string message)
{
	std::cout<< "\n"<< message << "\n";
}


class Model
{
private:

public:

};

//Class to solve problems in linear elasticity 
class ElasticSolids : public Model
{
private:
	double E, poisson, lambda, G; //Lame parameters.

	std::vector<std::vector<double>> global_k = {};

	Matrix elastic_matrix;
	Matrix b_matrix;

	void Lame()
	{
		lambda = (E*poisson) / ((1+poisson)*(1-(2*poisson)));
		G = E / (2 * (1 + poisson));
	}

	//Populates isotropic elasticity matrix with Lame parameters
	void ConstructElasticMatrix()
	{
		std::vector<std::vector<double>> temp = {
			{((2*G)+lambda),lambda,lambda,0,0,0},
			{lambda,((2 * G) + lambda),lambda,0,0,0},
			{lambda,lambda,((2 * G) + lambda),0,0,0},
			{0,0,0,G,0,0},
			{0,0,0,0,G,0},
			{0,0,0,0,0,G} };
		elastic_matrix = Matrix::Matrix(temp);
	}
	//Constructs 6x12 elemental B matrix
	//Requires global shape function derivatives
	void ConstructBMatrix(TetrahedralElement& el,Matrix& temp_b)
	{
		std::vector<std::vector<double>> nodes = el.GetNodes();
		std::vector<std::vector<double>> global_shape_derivatives = el.GetGlobalShapeDerivatives(0,0,0);

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
		temp_b = Matrix::Matrix(Mat);
	}

public:
	ElasticSolids()//Default constructors
	{
		std::cout << "Default constructor: ElasticSolids()" << std::endl;
	}

	ElasticSolids(Reader& reader, double e, double pois)
	{
		E = e;
		poisson = pois;
		Lame();
		ConstructElasticMatrix();

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

	Quadrature integrator;

	std::cout << integrator.Integrate(2,Funky);

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


}

