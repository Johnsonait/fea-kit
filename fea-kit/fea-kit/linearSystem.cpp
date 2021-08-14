#include <vector>
#include <cmath>
#include <iostream>

#include "linearsystem.h"


std::vector<std::vector<double>> A;
std::vector<double> b;

//Return length of vector, useful for calculating residuals
double LinearSystem::VectorNorm(std::vector<double>& x) 
{
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += std::pow(x[i], 2);
    }
    return std::pow(sum, 0.5);
}
//Return absolute value of a double, useful for calculating residuals
double LinearSystem::Abs(const double& val)
{
    return val < 0 ? -1 * val : val;
}

//Class constructor takes the square system matrix and constraint vector Ax = b --> Matrix*x = Vec
LinearSystem::LinearSystem(const std::vector<std::vector<double>>& Matrix, const std::vector<double>& Vec) 
    : A(Matrix), b(Vec) {}

//Solve system of equations, solution outputted to x that is passed in by reference
//Uses Gauss-siedel method
void LinearSystem::Solve(std::vector<double>& x)
{
    x.resize(b.size(), 0.1);

    double residual = 1.0;
    int count = 0;

    while (residual > 0.000001 && count < 10000)
    {
        double normStore = VectorNorm(x);
        for (int i = 0; i < b.size(); i++)
        {
            x[i] = b[i];
            for (int j = 0; j < b.size(); j++)
            {
                if (i != j)
                {
                    x[i] -= x[j] * A[i][j];
                }
            }
            x[i] *= 1 / A[i][i];
        }
        residual = abs(normStore - VectorNorm(x));
        count += 1;
    }
}

//Print solution vector (x) to console
void LinearSystem::PrintSol(std::vector<double>& x)
{
    for (int i = 0; i < x.size(); i++)
    {
        std::cout << "x" << i << ": " << x[i] << std::endl;
    }
}
