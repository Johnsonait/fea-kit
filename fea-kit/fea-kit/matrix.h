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
	Matrix operator * (const Matrix&); //Matrix multiplication
	Matrix operator * (const Matrix&) const; //Matrix multiplication for immutable object
	Matrix operator * (const double& a); //Scalar multiplication
	Matrix operator + (const Matrix&); //Matrix addition
	Matrix operator - (const Matrix&); //Matrix subtraction
	Matrix operator ^ (const size_t&); //Matrix powers (via repeated multiplication)
	std::vector<double>& operator [](const std::size_t&); //Data access operator
	const std::vector<double>& operator [](const std::size_t&) const; //Data access operator for immutable object


	//Useful
	Matrix& Transpose();
	Matrix GetTranspose() const;
	void PrintMatrix();
	void PrintMatrix() const;
	int CountRows();
	int CountRows() const;

	int CountCols();
	int CountCols() const;

	//Accessors 
	std::vector<std::vector<double>>& GetMatrix();
	const std::vector<std::vector<double>>& GetMatrix() const;
};
