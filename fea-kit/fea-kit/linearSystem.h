#pragma once
#include <vector>
#include <cmath>
#include <iostream>

#include "matrix.h"

class LinearSystem
{
private:

    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> b;

    double VectorNorm(std::vector<double>& x);
    double VectorNorm(std::vector<std::vector<double>>& x);

    double Abs(const double& val);

public:
    LinearSystem();
    LinearSystem(const std::vector<std::vector<double>>& mat,const std::vector<std::vector<double>>& vec);
    LinearSystem(std::shared_ptr<std::vector<std::vector<double>>>& mat, std::shared_ptr<std::vector<std::vector<double>>>& vec);
    LinearSystem(const Matrix& mat, const Matrix& vec);
    LinearSystem(Matrix& mat, Matrix& vec);

    void Solve(std::vector<double>& x);
    void Solve(std::shared_ptr<std::vector<std::vector<double>>> x);

    void PrintSol(std::vector<double>& x);
};
