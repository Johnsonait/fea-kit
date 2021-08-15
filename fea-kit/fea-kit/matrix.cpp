#include "matrix.h"


int Matrix::CountRows()
{
	return matrix.size();
}

int Matrix::CountRows() const 
{
	return matrix.size();
}

int Matrix::CountCols()
{
	return matrix[0].size();
}
int Matrix::CountCols() const
{
	return matrix[0].size();
}

//Constructors
Matrix::Matrix() //Default constructor
{
	matrix = {};
}
Matrix::Matrix(const size_t& rows, const size_t& cols)
{
	matrix.resize(rows);
	for (size_t i = 0; i < rows; i++)
	{
		matrix[i].resize(cols);
		for (size_t j = 0; j < cols; j++)
		{
			matrix[i][j] = 0;
		}
	}
}

Matrix::Matrix(std::vector<std::vector<double>>& mat) : matrix(mat) {}

//Operators
Matrix Matrix::operator * (const Matrix& B)
{
	std::vector<std::vector<double>> temp;
	temp.resize(this->CountRows());
	for (size_t n = 0; n < temp.size(); ++n)
	{
		temp[n].resize(B.CountCols());
	}

	for (size_t i = 0; i < this->CountRows(); i++)
	{
		for (size_t j = 0; j < B.CountCols(); j++)
		{
			double sum = 0;
			for (size_t k = 0; k < B.CountRows(); k++)
			{
				sum += matrix[i][k] * B[k][j];
			}
			temp[i][j] = sum;
		}
	}
	Matrix ret(temp);
	return ret;
}

//Operator to be used for scalar multiplication of matrices
Matrix Matrix::operator * (const double& a)
{
	std::vector<std::vector<double>> temp;
	for (size_t i = 0; i < this->CountRows(); ++i)
	{
		temp.push_back({});
		for (size_t j = 0; j < this->CountCols(); ++j)
		{
			temp[i].push_back(a * (*this)[i][j]);
		}
	}
	Matrix ret(temp);
	return ret;
}

Matrix Matrix::operator + (const Matrix& b)
{
	std::vector<std::vector<double>> temp;
	if (this->CountCols() == b.CountCols() && this->CountRows() == b.CountRows())
	{
		for (size_t i = 0; i < matrix.size(); ++i) //For each row i
		{
			temp.push_back({});
			for (size_t j = 0; j < matrix[0].size(); ++i) //For each column j
			{
				temp[i].push_back(matrix[i][j]+b[i][j]);
			}
		}
	}
	Matrix ret(temp);
	return ret;
}

Matrix Matrix::operator - (const Matrix& b)
{
	std::vector<std::vector<double>> temp;
	if (this->CountCols() == b.CountCols() && this->CountRows() == b.CountRows())
	{
		for (size_t i = 0; i < matrix.size(); ++i) //For each row i
		{
			temp.push_back({});
			for (size_t j = 0; j < matrix[0].size(); ++i) //For each column j
			{
				temp[i].push_back(matrix[i][j] - b[i][j]);
			}
		}
	}
	Matrix ret(temp);
	return ret;
}

Matrix Matrix::operator ^ (const size_t& n)
{
	Matrix temp = *this;
	for (size_t i = 0; i < n; ++i)
	{
		temp = temp * (*this);
	}
	return temp;
}

std::vector<double>& Matrix::operator [] (const std::size_t& n)
{
	return matrix[n];
}

const std::vector<double>& Matrix::operator [] (const std::size_t& n) const
{
	return matrix[n];
}

Matrix& Matrix::Transpose()
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
	return *this;
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

std::vector<std::vector<double>>& Matrix::GetMatrix()
{
	return matrix;
}

const std::vector<std::vector<double>>& Matrix::GetMatrix() const
{

	return matrix;
}
