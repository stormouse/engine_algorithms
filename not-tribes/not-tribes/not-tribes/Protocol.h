#pragma once
#ifndef __NOTTRIBES_PROTOCOL_H__
#define __NOTTRIBES_PROTOCOL_H__

#include <stdint.h>
#include <assert.h>
#include <stdio.h>

namespace NTribes
{

	class Stream
	{
	public:

		uint8_t* GetBuffer() { return data; }

		int SizeBytes() const { return (bitIndex + 7) >> 3; }
		
		int SizeBits() const { return bitIndex; }

		virtual ~Stream() {}

	protected:

		Stream(int sizeInBytes, uint8_t* buf = nullptr)
			: sizeBits(sizeInBytes * 8)
			, data(buf)
			, byteIndex(0)
			, bitIndex(0)
		{}

		uint8_t* data;
		int sizeBits;
		int byteIndex;
		int bitIndex;
		inline int sizeBytes() { return (bitIndex + 7) >> 3; }
		inline void AlignByteIndex()
		{
			if (byteIndex * 8 < bitIndex)
			{
				byteIndex++;
				bitIndex = byteIndex * 8;
			}
		}
	};


	class WriteStream : public Stream
	{
	public:
		enum { IsReading = 0 };
		enum { IsWriting = 1 };

		WriteStream(int totalBytes) : Stream(totalBytes, new uint8_t[totalBytes]) 
		{
			for (int i = 0; i < totalBytes; i++)
			{
				data[i] = uint8_t(0);
			}
		}

		void WriteBits(uint32_t value, int bits)
		{
			assert(bits <= 32);
			assert(bits > 0);
			assert(bitIndex + bits <= sizeBits);

			// little endian
			value &= (uint64_t(1) << bits) - 1;
			int bitLeft = 8 - (bitIndex - (byteIndex * 8));
			while (bitLeft < bits)
			{
				uint8_t res = value >> (bits - bitLeft);
				data[byteIndex] |= res;
				bitIndex += bitLeft;
				bits -= bitLeft;
				value &= (uint64_t(1) << bits) - 1;
				bitLeft = 8;
				byteIndex++;
			}

			data[byteIndex] |= value << (bitLeft - bits);
			bitIndex += bits;
		}

		void WriteInt32(uint32_t value, bool align)
		{
			if (align)
			{
				assert(sizeBytes() + 4 <= sizeBits * 8);
				AlignByteIndex();

				// little endian
				data[byteIndex++] = (value >> 24) & 0xFF;
				data[byteIndex++] = (value >> 16) & 0xFF;
				data[byteIndex++] = (value >> 8) & 0xFF;
				data[byteIndex++] = 0xFF;

				bitIndex += 32;
			}
			else
			{
				WriteBits(value, 32);
			}
		}

		void WriteInt64(uint64_t value, bool align)
		{
			if (align)
			{
				assert(sizeBytes() + 8 <= sizeBits * 8);
				AlignByteIndex();

				// little endian
				data[byteIndex++] = (value >> 56) & 0xFF;
				data[byteIndex++] = (value >> 48) & 0xFF;
				data[byteIndex++] = (value >> 40) & 0xFF;
				data[byteIndex++] = (value >> 32) & 0xFF;
				data[byteIndex++] = (value >> 24) & 0xFF;
				data[byteIndex++] = (value >> 16) & 0xFF;
				data[byteIndex++] = (value >> 8) & 0xFF;
				data[byteIndex++] = 0xFF;

				bitIndex += 64;
			}
			else
			{
				// little endian
				WriteBits((value >> 32) & 0xFFFF, 32);
				WriteBits(value & 0xFFFF, 32);
			}
		}

		void WriteByte(uint8_t value, bool align)
		{
			if (align)
			{
				assert(sizeBytes() + 1 <= sizeBits * 8);
				AlignByteIndex();

				data[byteIndex++] = value;
				bitIndex += 8;
			}
			else
			{
				WriteBits((uint32_t)value, 8);
			}
		}

		void WriteData(uint8_t* data, int length, bool align)
		{
			if (align)
			{
				assert(sizeBytes() + length <= sizeBits * 8);
				AlignByteIndex();
			}
			else
			{
				assert(bitIndex + length * 8 <= sizeBits);
			}

			for (int i = 0; i < length; i++)
			{
				WriteByte(data[i], align);
			}
		}

		~WriteStream()
		{
			delete [] data;
		}

	};


	class ReadStream : public Stream
	{

	public:

		enum { IsReading = 1 };
		enum { IsWriting = 0 };

		ReadStream(uint8_t* data, int totalBytes) : Stream(totalBytes, data) {}

		uint32_t ReadBits(int bits)
		{
			assert(bits < 32);
			assert(bits > 0);
			assert(bitIndex + bits <= sizeBits);

			uint32_t value = 0;

			// little endian
			int bitLeft = 8 - (bitIndex - (byteIndex * 8));
			int mask = (uint64_t(1) << bitLeft) - 1;
			while (bitLeft < bits)
			{
				value |= (data[byteIndex] & mask) << (bits - bitLeft);
				bitIndex += bitLeft;
				bits -= bitLeft;
				bitLeft = 8;
				mask = 0xFF;
				byteIndex++;
			}

			value |= (data[byteIndex] & mask) >> (bitLeft - bits);
			bitIndex += bits;

			return value;
		}

		int32_t ReadInt32(bool align)
		{
			int32_t value = 0;

			if (align)
			{
				assert(sizeBytes() + 4 <= sizeBits * 8);
				AlignByteIndex();

				// little endian
				value |= data[byteIndex++] << 24;
				value |= data[byteIndex++] << 16;
				value |= data[byteIndex++] << 8;
				value |= data[byteIndex++];

				bitIndex += 32;
			}
			else
			{
				value = ReadBits(32);
			}

			return value;
		}

		int64_t WriteInt64(bool align)
		{
			int64_t value = 0;

			if (align)
			{
				assert(sizeBytes() + 8 <= sizeBits * 8);
				AlignByteIndex();

				// little endian
				value |= (data[byteIndex++] >> 56) & 0xFF;
				value |= (data[byteIndex++] >> 48) & 0xFF;
				value |= (data[byteIndex++] >> 40) & 0xFF;
				value |= (data[byteIndex++] >> 32) & 0xFF;
				value |= (data[byteIndex++] >> 24) & 0xFF;
				value |= (data[byteIndex++] >> 16) & 0xFF;
				value |= (data[byteIndex++] >> 8) & 0xFF;
				value |= data[byteIndex++] & 0xFF;

				bitIndex += 64;
			}
			else
			{
				// little endian
				value |= uint64_t(ReadBits(32)) << 32;
				value |= ReadBits(32);
			}

			return value;
		}

		uint8_t ReadByte(bool align)
		{
			uint8_t value = 0;

			if (align)
			{
				assert(sizeBytes() + 1 <= sizeBits * 8);
				AlignByteIndex();

				value = data[byteIndex++];
				bitIndex += 8;
			}
			else
			{
				value = ReadBits(8);
			}

			return value;
		}

		void ReadData(uint8_t* data, int length, bool align)
		{
			if (align)
			{
				assert(sizeBytes() + length <= sizeBits * 8);
				AlignByteIndex();
			}
			else
			{
				assert(bitIndex + length * 8 <= sizeBits);
			}

			for (int i = 0; i < length; i++)
			{
				data[i] = ReadByte(align);
			}
		}

		~ReadStream() {}
	};
}


#endif