#pragma once
#include <vector>
#include <iostream>

class Matrix
{
private:
	std::vector<std::vector<double>> matrix;

public:
	//Constructors
	Matrix();
	Matrix(const size_t& rows, const size_t& cols);
	Matrix(std::vector<std::vector<double>>& mat);

	//Operators
	Matrix operator * (const Matrix&);
	Matrix operator * (const double& a);
	Matrix operator + (const Matrix&);
	Matrix operator - (const Matrix&);
	Matrix operator ^ (const size_t&);
	std::vector<double>& operator [](const std::size_t&);
	const std::vector<double>& operator [](const std::size_t&) const;


	//Useful
	Matrix& Transpose();
	void PrintMatrix();
	int CountRows();
	int CountRows() const;

	int CountCols();
	int CountCols() const;

	//Accessors 

	std::vector<std::vector<double>>& GetMatrix();
	const std::vector<std::vector<double>>& GetMatrix() const;
};
