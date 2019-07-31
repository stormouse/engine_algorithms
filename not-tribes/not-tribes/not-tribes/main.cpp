#include <iostream>

#include "Protocol.h"

int main()
{
	char data[] = "Hello, world";

	NTribes::WriteStream writeStream(1024);
	writeStream.WriteBits((int)0x41424344, 16);
	writeStream.WriteBits((int)0x43444546, 16);
	writeStream.WriteBits((int)0x45464748, 16);
	writeStream.WriteData((uint8_t*)data, 13, false);

	NTribes::ReadStream readStream(writeStream.GetBuffer(), writeStream.SizeBytes());

	char buf[1024];
	readStream.ReadData((uint8_t*)buf, 19, true);

	std::cout << buf << std::endl;

	return 0;
}