#include <vector>
#include <string>

#include "body.h"


Body::Body()
{

}

void Body::GetStrain()
{

}

void Body::GetStress()
{

}

void Body::OutputData()
{

}
void Body::AddNode(double x,double y, double z)
{
	nodes.push_back({x,y,z});
}
