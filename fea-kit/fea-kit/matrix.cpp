#include <vector>
#include <iostream>

#include "matrix.h"


int Matrix::CountRows()
{
	return matrix.size();
}

int Matrix::CountCols()
{
	return matrix[0].size();
}

Matrix::Matrix() //Default constructor
{
}

Matrix::Matrix(std::vector<std::vector<double>>& mat) : matrix(mat) {}

Matrix Matrix::operator * (Matrix& B)
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

void Matrix::Transpose()
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

void Matrix::PrintMatrix()
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