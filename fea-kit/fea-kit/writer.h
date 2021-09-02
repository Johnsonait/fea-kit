#pragma once
#include <iostream>
#include <fstream>

#include "body.h"

class Writer
{
private:
	std::string file_name;
	Body* body_ptr;
public:
	Writer();
	Writer(const std::string& file, Body& body);

	void Write();
};

