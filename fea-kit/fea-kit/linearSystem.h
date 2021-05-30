#pragma once


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include "linearSystem.h"

class linearSystem
{
private:

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    double vectorNorm(std::vector<double>& x);

    double abs(double x);

public:

    linearSystem(std::vector<std::vector<double>> Matrix, std::vector<double> Vec);

    void solve(std::vector<double>& x);

    void printSol(std::vector<double>& x);
};

