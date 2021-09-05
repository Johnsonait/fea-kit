#include "program.h"

Program::Program() = default;

void Program::Run()
{
	Body body;

	Reader reader("model.txt",body);

	LinearElasticSolids model(reader,body);
	model.Solve();
	model.PostProcess();
	
	Writer writer("C:\\Users\\johan\\Desktop\\FEM\\fea-kit\\render\\results.txt", body);
	writer.Write();
}
