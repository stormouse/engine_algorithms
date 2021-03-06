// udp-test-program.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include "../reliable-udp/udp_network.h"


int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cout << "need argument --server or --client\n";
		return -1;
	}

	if (strcmp(argv[1], "--server") == 0) {
		return unet::ReceiverProgram();
	}

	if (strcmp(argv[1], "--client") == 0) {
		return unet::SenderProgram();
	}

	return 0;
}
