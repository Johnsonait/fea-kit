#include "linearsystem.h"

static const int COUNTMAX = 10000;
static const double MINRES = 1e-5;

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
double LinearSystem::VectorNorm(std::vector<std::vector<double>>& x)
{
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += std::pow(x[i][0], 2);
    }
    return std::pow(sum, 0.5);
}
//Return absolute value of a double, useful for calculating residuals
double LinearSystem::Abs(const double& val)
{
    return val < 0 ? -1 * val : val;
}

LinearSystem::LinearSystem() 
{
    matrix = { {} };
    b = { {} };
}

//Class constructor takes the square system matrix and constraint vector Ax = b --> Matrix*x = Vec
LinearSystem::LinearSystem(const std::vector<std::vector<double>>& mat, const std::vector<std::vector<double>>& vec)
{
    
    matrix = mat;
    b = vec;

}

LinearSystem::LinearSystem(std::shared_ptr<std::vector<std::vector<double>>>& mat, std::shared_ptr<std::vector<std::vector<double>>>& vec)
{
    matrix = *mat;
    b = *vec;
}

LinearSystem::LinearSystem(const Matrix& mat, const Matrix& vec)
{
    matrix = mat.GetMatrix();
    b = vec.GetMatrix();
}

LinearSystem::LinearSystem(Matrix& mat, Matrix& vec)
{
    matrix = mat.GetMatrix();
    b = vec.GetMatrix();
}

//Solve system of equations, solution outputted to x that is passed in by reference
//Uses Gauss-siedel method
void LinearSystem::Solve(std::vector<double>& x)
{
    x.resize(b[0].size(), 0.1);

    double residual = 1.0;
    int count = 0;

    while (residual > 0.000001 && count < 10000)
    {
        double normStore = VectorNorm(x);
        for (int i = 0; i < b[0].size(); i++)
        {
            x[i] = b[0][i];
            for (int j = 0; j < b[0].size(); j++)
            {
                if (i != j)
                {
                    x[i] -= x[j] * matrix[i][j];
                }
            }
            x[i] *= 1 / matrix[i][i];
        }
        residual = abs(normStore - VectorNorm(x));
        count += 1;
    }
}
void LinearSystem::Solve(std::shared_ptr<std::vector<std::vector<double>>> x)
{
    (*x)[0].resize(b[0].size(), 0.2);

    double residual = 1.0;
    double residual_store = 0;
    int count = 0;

    while (residual > MINRES && residual != residual_store && count < COUNTMAX)
    {
        double normStore = VectorNorm((*x)[0]);
        residual_store = residual;
        if (count % 1000 == 0)
        {
            std::cout << "Iterations: " << count << "\n";
            std::cout << "Residual: " << residual << "\n";
        }
        for (int i = 0; i < b.size(); i++)
        {
            (*x)[i][0] = b[i][0];
            for (int j = 0; j < b.size(); j++)
            {
                if (i != j)
                {
                    (*x)[i][0] -= (*x)[j][0] * matrix[i][j];
                }
            }
            (*x)[i][0] = (*x)[i][0] / matrix[i][i];
        }
        residual = abs(normStore - VectorNorm((*x)));
        count += 1;
    }
    if (count == COUNTMAX)
    {
        std::cout<< std::endl << "Max iterations reached, solution did not converge! Residual: "<< residual << std::endl;
        return;
    }
    std::cout << "Maximum Iterations: " << COUNTMAX << std::endl;
    std::cout << "Iterations: " << count << std::endl;
    std::cout << "Minimum Residual: " << MINRES << std::endl;
    std::cout << "Residual: " << residual << std::endl;
}

//Print solution vector (x) to console
void LinearSystem::PrintSol(std::vector<double>& x)
{
    for (int i = 0; i < x.size(); i++)
    {
        std::cout << "x" << i << ": " << x[i] << std::endl;
    }
}
