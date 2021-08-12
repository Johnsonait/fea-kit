#pragma once

class Matrix
{
private:
	std::vector<std::vector<double>> matrix;

	int CountRows();

	int CountCols();

public:
	//Constructors
	Matrix();

	Matrix(std::vector<std::vector<double>>& mat);

	//Operators
	Matrix operator * (Matrix& B);

	//Useful
	void Transpose();
	void PrintMatrix();

	//Accessors 

	std::vector<std::vector<double>>& GetMatrix();
};