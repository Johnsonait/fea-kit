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
#include "tetrahedral_element.h"
#include "program.h"
#include "model.h"
#include "linear_elastic_solids.h"

void Log(const std::string& message)
{
	std::cout<< "\n"<< message << "\n";
}

//Quick function to test Quadrature class
double Funky(double x, double y, double z)
{
	return x * x * y * y * z * z;
}

int main()
{
	/**/
	std::vector<std::vector<double>> temp_Mat = { {1,0,0,0},{2,4,3,0.5},{1,3,5,2},{1,1,1,1} };
	std::vector<std::vector<double>> temp_b = {{ 1,1,1,1 }};

	std::vector<double> x;
	 
	Matrix Mat(temp_Mat), b(temp_b);

	LinearSystem system(Mat, b);

	system.Solve(x);
	system.PrintSol(x);
	

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

	Matrix D = C * 1.5;
	D.PrintMatrix();
	D.Transpose();
	D.PrintMatrix();

	std::cout << "C: " << std::endl;
	C.PrintMatrix();
	std::cout << "D: " << std::endl;
	D.PrintMatrix();

	(D*C).PrintMatrix();
	

	LinearElasticSolids model;
	model.Solve();
}

