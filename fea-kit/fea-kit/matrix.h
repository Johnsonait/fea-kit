#pragma once

class Matrix
{
private:
	std::vector<std::vector<double>> matrix;

	int CountRows();

	int CountCols();

public:
	Matrix();

	Matrix(std::vector<std::vector<double>>& mat);

	Matrix operator * (Matrix& B);

	void Transpose();

	void PrintMatrix();
};