#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include "mesh.h"
#include "linearSystem.h"


int main()
{
    mesh MESH("C:\\Users\\johan\\Desktop\\FEM\\fea-kit\\mesh.txt");

    //MESH.printNodes();
    //MESH.printElements();
    //MESH.printShellList();
    std::vector<std::vector<double>> Mat = { {1,0,0},{0,1,0},{0,0,1} };
    std::vector<double> Con = { 1,-2,1 };
    std::vector<double> Sol;

    linearSystem testSystem(Mat, Con);
    testSystem.solve(Sol);
    testSystem.printSol(Sol);

    return 0;
}