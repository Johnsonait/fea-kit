#pragma once

class LinearSystem
{
private:

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    double VectorNorm(std::vector<double>& x);

    double Abs(const double& val);

public:

    LinearSystem(std::vector<std::vector<double>>& Matrix, std::vector<double>& Vec);

    void Solve(std::vector<double>& x);

    void PrintSol(std::vector<double>& x);
};