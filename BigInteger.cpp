#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <chrono>
#include <string>
#include<stdexcept>
#include<cstring>
#include <memory>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#include <vector>
#include <memory>

using namespace std;

class IIRandom;

typedef shared_ptr<IIRandom> IRandom;

class IIRandom
{
public:
	virtual void NextBytes(vector<uint8_t> &buf) = 0;

	virtual double NextDouble() = 0;

	virtual int32_t Next() = 0;
	virtual int32_t Next(const int32_t maxValue) = 0;
	virtual int32_t Next(const int32_t minValue, const int32_t maxValue) = 0;

};

#include <vector>

using namespace std;

class IIRandomGenerator;

typedef shared_ptr<IIRandomGenerator> IRandomGenerator;

class IIRandomGenerator
{
public:
	/// <summary>Add more seed material to the generator.</summary>
	/// <param name="seed">A byte array to be mixed into the generator's state.</param>
	virtual void AddSeedMaterial(vector<uint8_t> &seed) = 0;

	/// <summary>Add more seed material to the generator.</summary>
	/// <param name="seed">A long value to be mixed into the generator's state.</param>
	virtual void AddSeedMaterial(const int64_t seed) = 0;

	/// <summary>Fill byte array with random values.</summary>
	/// <param name="bytes">Array to be filled.</param>
	virtual void NextBytes(vector<uint8_t> &bytes) = 0;

	/// <summary>Fill byte array with random values.</summary>
	/// <param name="bytes">Array to receive bytes.</param>
	/// <param name="start">Index to start filling at.</param>
	/// <param name="len">Length of segment to fill.</param>
	virtual void NextBytes(vector<uint8_t> &bytes, const int32_t start, const int32_t len) = 0;
};

#include <vector>
#include <memory>

using namespace std;

class IIRandomNumberGenerator;

typedef shared_ptr<IIRandomNumberGenerator> IRandomNumberGenerator;

class IIRandomNumberGenerator
{
public:
	virtual void GetBytes(vector<uint8_t> &data) = 0;

	virtual void GetNonZeroBytes(vector<uint8_t> &data) = 0;

};

class IICryptoApiRandomGenerator;

typedef shared_ptr<IICryptoApiRandomGenerator> ICryptoApiRandomGenerator;

class IICryptoApiRandomGenerator : public virtual IIRandomGenerator
{};

class IISecureRandom;

typedef shared_ptr<IISecureRandom> ISecureRandom;

class IISecureRandom : public virtual IIRandom
{
public:
	virtual vector<uint8_t> GenerateSeed(const int32_t length) = 0;
	virtual void SetSeed(vector<uint8_t> &seed) = 0;
	virtual void SetSeed(const int64_t seed) = 0;

	virtual void NextBytes(vector<uint8_t> &buf) = 0;
	virtual void NextBytes(vector<uint8_t> &buf, const int32_t off, const int32_t len) = 0;
	virtual int32_t NextInt32() = 0;
	virtual int64_t NextInt64() = 0;

};

class NumberStyles
{
public:
	static const auto None = 0;
	static const auto AllowLeadingWhite = 1;
	static const auto AllowTrailingWhite = 2;
	static const auto AllowLeadingSign = 4;
	static const auto AllowTrailingSign = 8;
	static const auto AllowParentheses = 16;
	static const auto AllowDecimalPoint = 32;
	static const auto AllowThousands = 64;
	static const auto AllowExponent = 128;
	static const auto AllowCurrencySymbol = 256;
	static const auto AllowHexSpecifier = 512;
	static const auto Integer = 4 | 2 | 1;

};

class Bits
{
public:
	static void ReverseByteArray(const void *Source, void * Dest, int64_t size)
	{
		uint8_t *ptr_src = (uint8_t *)Source;
		uint8_t *ptr_dest = (uint8_t *)Dest;

		ptr_dest = ptr_dest + (size - 1);
		while (size > 0)
		{
			*ptr_dest = *ptr_src;
			ptr_src += 1;
			ptr_dest -= 1;
			size -= 1;
		} // end while
	} // end function ReverseByteArray

	inline static int32_t ReverseBytesInt32(const int32_t value)
	{
		int32_t i1 = value & 0xFF;
		int32_t i2 = Bits::Asr32(value, 8) & 0xFF;
		int32_t i3 = Bits::Asr32(value, 16) & 0xFF;
		int32_t i4 = Bits::Asr32(value, 24) & 0xFF;

		return (i1 << 24) | (i2 << 16) | (i3 << 8) | (i4 << 0);
	} // end function ReverseBytesInt32

	inline static uint8_t ReverseBitsUInt8(const uint8_t value)
	{
		uint8_t result = ((value >> 1) & 0x55) | ((value << 1) & 0xAA);
		result = ((result >> 2) & 0x33) | ((result << 2) & 0xCC);
		return ((result >> 4) & 0x0F) | ((result << 4) & 0xF0);
	} // end function ReverseBitsUInt8

	inline static uint16_t ReverseBytesUInt16(const uint16_t value)
	{
		return ((value & uint32_t(0xFF)) << 8 | (value & uint32_t(0xFF00)) >> 8);
	} // end function ReverseBytesUInt16

	inline static uint32_t ReverseBytesUInt32(const uint32_t value)
	{
		return (value & uint32_t(0x000000FF)) << 24 |
			(value & uint32_t(0x0000FF00)) << 8 |
			(value & uint32_t(0x00FF0000)) >> 8 |
			(value & uint32_t(0xFF000000)) >> 24;
	} // end function ReverseBytesUInt32

	inline static uint64_t ReverseBytesUInt64(const uint64_t value)
	{
		return (value & uint64_t(0x00000000000000FF)) << 56 |
			(value & uint64_t(0x000000000000FF00)) << 40 |
			(value & uint64_t(0x0000000000FF0000)) << 24 |
			(value & uint64_t(0x00000000FF000000)) << 8 |
			(value & uint64_t(0x000000FF00000000)) >> 8 |
			(value & uint64_t(0x0000FF0000000000)) >> 24 |
			(value & uint64_t(0x00FF000000000000)) >> 40 |
			(value & uint64_t(0xFF00000000000000)) >> 56;
	} // end function ReverseBytesUInt64

	inline static int32_t Asr32(const int32_t value, const int32_t ShiftBits)
	{
		int32_t result = value >> ShiftBits;
		if ((value & 0x80000000) > 0)
		{
			// if you don't want to cast (0xFFFFFFFF) to an Int32,
			// simply replace it with (-1) to avoid range check error.
			result = result | (int32_t(0xFFFFFFFF) << (32 - ShiftBits));
		} // end if

		return result;
	} // end function Asr32

	inline static int64_t Asr64(const int64_t value, const int32_t ShiftBits)
	{
		int64_t result = value >> ShiftBits;
		if ((value & 0x8000000000000000) > 0)
		{
			result = result | (0xFFFFFFFFFFFFFFFF << (64 - ShiftBits));
		} // end if

		return result;
	} // end function Asr64

	inline static uint32_t RotateLeft32(const uint32_t a_value, int32_t a_n)
	{
		a_n = a_n & 31;
		return (a_value << a_n) | (a_value >> (32 - a_n));
	} // end function RotateLeft32

	inline static uint64_t RotateLeft64(const uint64_t a_value, int32_t a_n)
	{
		a_n = a_n & 63;
		return (a_value << a_n) | (a_value >> (64 - a_n));
	} // end function RotateLeft64

	inline static uint32_t RotateRight32(const uint32_t a_value, int32_t a_n)
	{
		a_n = a_n & 31;
		return (a_value >> a_n) | (a_value << (32 - a_n));
	} // end function RotateRight32

	inline static uint64_t RotateRight64(const uint64_t a_value, int32_t a_n)
	{
		a_n = a_n & 63;
		return (a_value >> a_n) | (a_value << (64 - a_n));
	} // end function RotateRight64


};

#include <memory>
#include <chrono>

using namespace std::chrono;

class Pcg
{
private:
	/// <summary>
	/// The RNG state. All values are possible.
	/// </summary>
	static uint64_t m_state;

	/// <summary>
	/// Controls which RNG sequence (stream) is selected.
	/// Must <strong>always</strong> be odd.
	/// </summary>
	static uint64_t m_inc;

	/// <summary>
	/// Internal variable used for Casting.
	/// </summary>
	static int64_t I64;

	/// <summary>
	/// Seed Pcg in two parts, a state initializer
	/// and a sequence selection constant (a.k.a.
	/// stream id).
	/// </summary>
	/// <param name="initState">Initial state.</param>
	/// <param name="initSeq">Initial sequence</param>
	static inline void Seed(const uint64_t initState, const uint64_t initSeq);

	/// <summary>
	/// Generates a uniformly distributed number, r,
	/// where 0 &lt;= r &lt; exclusiveBound.
	/// </summary>
	/// <param name="exclusiveBound">Exclusive bound.</param>
	static uint32_t Range32(const uint32_t exclusiveBound);

	/// <summary>
	/// Generates an Init State from System Time.
	/// </summary>
	/// <param name="initSeq">Calculated initSeq.</param>
	/// <returns> UInt64 </returns>
	static inline uint64_t GetInitState(uint64_t &initSeq);

	/// <summary>
	/// Generates an Init Sequence from GetInitState value * 181783497276652981.
	/// <param name="tempVal">Previous value from GetInitState.</param>
	/// </summary>
	/// <returns> UInt64 </returns>
	static inline uint64_t GetInitSeq(const uint64_t tempVal);

public:
	/// <summary>
	/// Initializes a new instance of the <see cref="TPcg"/> class
	/// <strong>FOR TESTING</strong> with a <strong>KNOWN</strong> seed.
	/// </summary>
	Pcg();

	/// <summary>
	/// Initializes a new instance of the <see cref="TPcg"/> class.
	/// </summary>
	/// <param name="initState">Initial state.</param>
	/// <param name="initSeq">Initial sequence</param>
	Pcg(const uint64_t initState, const uint64_t initSeq);

	/// <summary>
	/// Generates a uniformly-distributed 32-bit random number.
	/// </summary>
	static uint32_t NextUInt32();

	/// <summary>
	/// Generates a uniformly distributed number, r,
	/// where minimum &lt;= r &lt; exclusiveBound.
	/// </summary>
	/// <param name="minimum">The minimum inclusive value.</param>
	/// <param name="exclusiveBound">The maximum exclusive bound.</param>
	static uint32_t NextUInt32(const uint32_t minimum, const uint32_t exclusiveBound);

	/// <summary>
	/// Generates a uniformly distributed number, r,
	/// where minimum &lt;= r &lt; exclusiveBound.
	/// </summary>
	/// <param name="minimum">The minimum inclusive value.</param>
	/// <param name="exclusiveBound">The maximum exclusive bound.</param>
	static int NextInt(const int minimum, const int exclusiveBound);

};

uint64_t Pcg::m_state = 0;
uint64_t Pcg::m_inc = 0;
int64_t Pcg::I64 = 0;

Pcg::Pcg()
{
	uint64_t LinitState, LinitSeq;

	LinitState = GetInitState(LinitSeq);

	// ==> initializes using system time as initState and calculated value as
	// initSeq
	Seed(LinitState, LinitSeq);
}

Pcg::Pcg(const uint64_t initState, const uint64_t initSeq)
{
	Seed(initState, initSeq);
}

uint32_t Pcg::NextUInt32()
{
	uint64_t oldState;
	uint32_t xorShifted;
	int rot;

	oldState = m_state;
	m_state = oldState * uint64_t(6364136223846793005) + m_inc;
	xorShifted = uint32_t(((oldState >> 18) ^ oldState) >> 27);
	rot = int(oldState >> 59);

	return (xorShifted >> rot) | (xorShifted << ((-rot) & 31));
}

uint32_t Pcg::NextUInt32(const uint32_t minimum, const uint32_t exclusiveBound)
{
	uint32_t boundRange, rangeResult;

	boundRange = exclusiveBound - minimum;
	rangeResult = Range32(boundRange);

	return rangeResult + minimum;
}

int Pcg::NextInt(const int minimum, const int exclusiveBound)
{
	uint32_t boundRange, rangeResult;

	boundRange = uint32_t(exclusiveBound - minimum);
	rangeResult = Range32(boundRange);

	return int(rangeResult) + int(minimum);
}

void Pcg::Seed(const uint64_t initState, const uint64_t initSeq)
{
	m_state = uint32_t(0);
	m_inc = (initSeq << 1) | uint64_t(1);
	NextUInt32();
	m_state = m_state + initState;
	NextUInt32();
}

uint32_t Pcg::Range32(const uint32_t exclusiveBound)
{
	uint32_t r, threshold;

	// To avoid bias, we need to make the range of the RNG
	// a multiple of bound, which we do by dropping output
	// less than a threshold. A naive scheme to calculate the
	// threshold would be to do
	//
	// threshold = UInt64($100000000) mod exclusiveBound;
	//
	// but 64-bit div/mod is slower than 32-bit div/mod
	// (especially on 32-bit platforms). In essence, we do
	//
	// threshold := UInt32((UInt64($100000000) - exclusiveBound) mod exclusiveBound);
	//
	// because this version will calculate the same modulus,
	// but the LHS value is less than 2^32.
	threshold = uint32_t((uint64_t(0x100000000) - exclusiveBound) % exclusiveBound);

	// Uniformity guarantees that this loop will terminate.
	// In practice, it should terminate quickly; on average
	// (assuming all bounds are equally likely), 82.25% of
	// the time, we can expect it to require just one
	// iteration. In the worst case, someone passes a bound
	// of 2^31 + 1 (i.e., 2147483649), which invalidates
	// almost 50% of the range. In practice bounds are
	// typically small and only a tiny amount of the range
	// is eliminated.
	while (true)
	{
		r = NextUInt32();
		if (r >= threshold)
		{
			return r % exclusiveBound;
		}
	}

	return  0; // to make FixInsight Happy :)
}

uint64_t Pcg::GetInitState(uint64_t &initSeq)
{
	milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());

	uint64_t result = ms.count();
	initSeq = GetInitSeq(result) * uint64_t(int64_t(1000000));

	return result;
}

uint64_t Pcg::GetInitSeq(const uint64_t tempVal)
{
	return tempVal * uint64_t(181783497276652981);
}

#include <string>

class RandomNumberGenerator : public IRandomNumberGenerator
{
private:
	// Resource string
	static const char * UnknownAlgorithm;

public:
	static IRandomNumberGenerator CreateRNG();

	static IRandomNumberGenerator CreateRNG(const string rngName);

	virtual void GetBytes(vector<uint8_t> &data) = 0;

	virtual void GetNonZeroBytes(vector<uint8_t> &data) = 0;

};

class PCGRandomNumberGenerator : public virtual IIRandomNumberGenerator, public RandomNumberGenerator
{
public:
	PCGRandomNumberGenerator();

	virtual void GetBytes(vector<uint8_t> &data);

	virtual void GetNonZeroBytes(vector<uint8_t> &data);

};

PCGRandomNumberGenerator::PCGRandomNumberGenerator()
{
	Pcg();
}

void PCGRandomNumberGenerator::GetBytes(vector<uint8_t> &data)
{
	int64_t i;

	for (i = data.size() - 1; i >= 0; i--)
	{
		data[i] = uint8_t(Pcg::NextInt(INT32_MIN, INT32_MAX));
	}
}

void PCGRandomNumberGenerator::GetNonZeroBytes(vector<uint8_t> &data)
{
	int64_t i;
	uint8_t val;

	i = data.size();
	while (i > 0)
	{
		do
		{
			val = uint8_t(Pcg::NextUInt32(INT32_MIN, INT32_MAX));
		} while (!(val == 0));

		data[i - 1] = val;
		i--;
	}
}

#include<stdexcept>
#include<cstring>

const char * RandomNumberGenerator::UnknownAlgorithm = "Unknown Random Generation Algorithm Requested";

IRandomNumberGenerator RandomNumberGenerator::CreateRNG()
{
	return make_shared<PCGRandomNumberGenerator>();
}

IRandomNumberGenerator RandomNumberGenerator::CreateRNG(const string rngName)
{
	if (rngName == "PCGRandomNumberGenerator") return make_shared<PCGRandomNumberGenerator>();

	throw invalid_argument(UnknownAlgorithm);
}

class CryptoApiRandomGenerator : public virtual IICryptoApiRandomGenerator, public virtual IIRandomGenerator
{
private:
	// Resource string
	const char * NegativeOffset = "Start Offset Cannot be Negative, \"Start\"";
	const char * ArrayTooSmall = "Byte Array Too Small For Requested Offset and Length";

	IRandomNumberGenerator FrndProv = nullptr;

public:
	/// <summary>
	/// Uses TRandomNumberGenerator.CreateRNG() to Get randomness generator
	/// </summary>
	CryptoApiRandomGenerator()
	{
		FrndProv = RandomNumberGenerator::CreateRNG();
	}

	CryptoApiRandomGenerator(IRandomNumberGenerator rng)
	{
		FrndProv = rng;
	}

	/// <summary>Add more seed material to the generator.</summary>
	/// <param name="seed">A byte array to be mixed into the generator's state.</param>
	virtual void AddSeedMaterial(vector<uint8_t> &seed) {};

	/// <summary>Add more seed material to the generator.</summary>
	/// <param name="seed">A long value to be mixed into the generator's state.</param>
	virtual void AddSeedMaterial(const int64_t seed) {};

	/// <summary>Fill byte array with random values.</summary>
	/// <param name="bytes">Array to be filled.</param>
	virtual void NextBytes(vector<uint8_t> &bytes)
	{
		FrndProv.get()->GetBytes(bytes);
	}

	/// <summary>Fill byte array with random values.</summary>
	/// <param name="bytes">Array to receive bytes.</param>
	/// <param name="start">Index to start filling at.</param>
	/// <param name="len">Length of segment to fill.</param>
	virtual void NextBytes(vector<uint8_t> &bytes, const int32_t start, const int32_t len)
	{
		vector<uint8_t> tmpBuf;

		if (start < 0)
			throw invalid_argument(NegativeOffset);

		if (bytes.size() < (start + len))
			throw invalid_argument(ArrayTooSmall);

		if ((bytes.size() == len) && (start == 0))
			NextBytes(bytes);
		else
		{
			tmpBuf = vector<uint8_t>(len);
			NextBytes(tmpBuf);

			std::memmove(&tmpBuf[0], &bytes[start], len * sizeof(uint8_t));
		}
	}

};

#include <memory>
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

class Random : public virtual IIRandom
{
private:
	// Resourcestring
	const char * BufferNil = "Buffer Cannot be Nil";
	const char * MaxValueNegative = "maxValue Must be Positive";
	const char * InvalidMinValue = "minValue Cannot be Greater Than maxValue";

private:
	const int32_t MBIG = int32_t(2147483647);
	const int32_t MSEED = int32_t(161803398);
	const int32_t MZ = int32_t(0);

	static vector<int32_t> SeedArray;

private:
	int32_t inext, inextp;

	inline int32_t InternalSample();
	double GetSampleForLargeRange();

protected:
	/// <summary>Returns a random floating-point number between 0.0 and 1.0.</summary>
	/// <returns>A double-precision floating point number that is greater than or equal to 0.0, and less than 1.0.</returns>
	virtual double Sample();

public:
	Random()
	{
		SeedArray = vector<int32_t>(56);
	}

	Random(const int32_t Seed);
	~Random();

	virtual void NextBytes(vector<uint8_t> &buf);

	virtual double NextDouble();

	virtual int32_t Next();
	virtual int32_t Next(const int32_t maxValue);
	virtual int32_t Next(const int32_t minValue, const int32_t maxValue);

};

vector<int32_t> Random::SeedArray = vector<int32_t>();

Random::Random(const int32_t Seed)
{
	Random();

	int32_t num1, num2, index1, index2;

	num1 = MSEED - abs(Seed);
	SeedArray[55] = num1;
	num2 = 1;
	for (index1 = 1; index1 < 55; index1++)
	{
		index2 = 21 * index1 % 55;
		SeedArray[index2] = num2;
		num2 = num1 - num2;
		if (num2 < 0)
			num2 = num2 + INT32_MAX;
		num1 = SeedArray[index2];
	}

	index1 = 1;
	while (index1 < 5)
	{
		for (index2 = 1; index2 < 56; index2++)
		{
			SeedArray[index2] = SeedArray[index2] - SeedArray[1 + (index2 + 30) % 55];
			if (SeedArray[index2] < 0)
				SeedArray[index2] = SeedArray[index2] + INT32_MAX;
		}

		index1++;
	}

	inext = 0;
	inextp = 21;
}

Random::~Random()
{
}

double Random::Sample()
{
	return InternalSample() * 4.6566128752458E-10;
}

void Random::NextBytes(vector<uint8_t> &buf)
{
	int32_t i;

	if (buf.empty())
		throw invalid_argument(BufferNil);

	for (i = 0; i < buf.size(); i++)
		buf[i] = uint8_t(InternalSample() % (255 + 1));
}

double Random::NextDouble()
{
	return Sample();
}

int32_t Random::InternalSample()
{
	int32_t _inext, _inextp, index1, index2, num;

	_inext = inext;
	_inextp = inextp;
	index1 = _inext + 1;
	if ((index1) >= 56)
		index1 = 1;

	index2 = _inextp + 1;
	if ((index2) >= 56)
		index2 = 1;

	num = SeedArray[index1] - SeedArray[index2];
	if (num < 0)
		num = num + INT32_MAX;

	SeedArray[index1] = num;
	inext = index1;
	inextp = index2;

	return num;
}

double Random::GetSampleForLargeRange()
{
	int32_t num;

	num = InternalSample();
	if (InternalSample() % 2 == 0)
		num = -num;

	return (num + 2147483646.0) / 4294967293.0;
}

int32_t Random::Next()
{
	return InternalSample();
}

int32_t Random::Next(const int32_t maxValue)
{
	if (maxValue < 0)
		throw out_of_range(MaxValueNegative);

	return int32_t(trunc(Sample() * maxValue));
}

int32_t Random::Next(const int32_t minValue, const int32_t maxValue)
{
	int64_t num;

	if (minValue > maxValue)
		throw out_of_range(InvalidMinValue);

	num = int64_t(maxValue) - int64_t(minValue);
	if (num <= int64_t(INT32_MAX))
		return int32_t(trunc(Sample()) * num) + minValue;

	return int32_t(int64_t(trunc(GetSampleForLargeRange()) * num) + int64_t(minValue));
}

class SecureRandom : public Random, public virtual IISecureRandom
{
private:
	// Resource string
	const char * UnrecognisedPRNGAlgorithm = "Unrecognised PRNG Algorithm: %s \"algorithm\"";
	const char * CannotBeNegative = "Cannot be Negative \"maxValue\"";
	const char * InvalidMaxValue = "maxValue Cannot be Less Than minValue";

private:
	static int64_t Counter;
	static ISecureRandom master;
	static double DoubleScale;

	static int64_t NextCounterValue();

protected:
	static IRandomGenerator generator;

public:
	SecureRandom()
	{
		Boot();
	}

	SecureRandom(IRandomGenerator _generator)
		: Random(0)
	{
		generator = _generator;
	}

	static ISecureRandom GetMaster()
	{
		return master;
	}

	virtual vector<uint8_t> GenerateSeed(const int32_t length);
	virtual void SetSeed(vector<uint8_t> &seed);
	virtual void SetSeed(const int64_t seed);

	virtual void NextBytes(vector<uint8_t> &buf);
	virtual void NextBytes(vector<uint8_t> &buf, const int32_t off, const int32_t len);
	virtual int32_t NextInt32();
	virtual int64_t NextInt64();

	virtual double NextDouble();

	virtual int32_t Next();
	virtual int32_t Next(const int32_t maxValue);
	virtual int32_t Next(const int32_t minValue, const int32_t maxValue);

	static vector<uint8_t> GetNextBytes(ISecureRandom secureRandom, const int32_t length);

	static void Boot();

};

#include <chrono>

using namespace std::chrono;

int64_t SecureRandom::Counter = 0;
double SecureRandom::DoubleScale = 0;
IRandomGenerator SecureRandom::generator = nullptr;
ISecureRandom SecureRandom::master = nullptr;


vector<uint8_t> SecureRandom::GenerateSeed(const int32_t length)
{
	return GetNextBytes(master, length);
}

void SecureRandom::SetSeed(vector<uint8_t> &seed)
{
	generator.get()->AddSeedMaterial(seed);
}

void SecureRandom::SetSeed(const int64_t seed)
{
	generator.get()->AddSeedMaterial(seed);
}

void SecureRandom::NextBytes(vector<uint8_t> &buf)
{
	generator.get()->NextBytes(buf);
}

void SecureRandom::NextBytes(vector<uint8_t> &buf, const int32_t off, const int32_t len)
{
	generator.get()->NextBytes(buf, off, len);
}

int32_t SecureRandom::NextInt32()
{
	uint32_t tempRes;
	vector<uint8_t> bytes(4);

	NextBytes(bytes);

	tempRes = bytes[0];
	tempRes = tempRes << 8;
	tempRes = tempRes | bytes[1];
	tempRes = tempRes << 8;
	tempRes = tempRes | bytes[2];
	tempRes = tempRes << 8;
	tempRes = tempRes | bytes[3];

	return int32_t(tempRes);
}

int64_t SecureRandom::NextInt64()
{
	return (int64_t(uint32_t(NextInt32())) << 32) | (int64_t(uint32_t(NextInt32())));
}

double SecureRandom::NextDouble()
{
	return uint64_t(NextInt64()) / DoubleScale;
}

int64_t SecureRandom::NextCounterValue()
{
	uint32_t LCounter;

	LCounter = uint32_t(Counter);
	int64_t result = LCounter; // InterLockedIncrement(LCounter);
	Counter = int64_t(LCounter);

	return result;
}

int32_t SecureRandom::Next()
{
	return NextInt32() & INT32_MAX;
}

int32_t SecureRandom::Next(const int32_t maxValue)
{
	int32_t bits;

	if (maxValue < 2)
	{
		if (maxValue < 0)
			throw out_of_range(CannotBeNegative);

		return 0;
	}


	// Test whether maxValue is a power of 2
	if ((maxValue & (maxValue - 1)) == 0)
	{
		bits = NextInt32() & INT32_MAX;
		return int32_t(Bits::Asr64((int64_t(bits) * maxValue), 31));
	}

	int32_t result;

	do
	{
		bits = NextInt32() & INT32_MAX;
		result = bits % maxValue;
		// Ignore results near overflow
	} while (((bits - (result + (maxValue - 1))) < 0));

	return result;
}

int32_t SecureRandom::Next(const int32_t minValue, const int32_t maxValue)
{
	int32_t diff, i;

	if (maxValue <= minValue)
	{
		if (maxValue == minValue) return minValue;

		throw invalid_argument(InvalidMaxValue);
	}


	diff = maxValue - minValue;
	if (diff > 0) return  minValue + Next(diff);

	while (true)
	{
		i = NextInt32();

		if ((i >= minValue) && (i < maxValue)) return i;
	}

	return 0; // to make FixInsight Happy :)
}

vector<uint8_t> SecureRandom::GetNextBytes(ISecureRandom secureRandom, const int32_t length)
{
	vector<uint8_t> result(length);
	secureRandom.get()->NextBytes(result);

	return result;
}

void SecureRandom::Boot()
{
	auto now = high_resolution_clock::now();
	nanoseconds nn = duration_cast<nanoseconds>(now.time_since_epoch());

	Counter = nn.count();
	if (master == nullptr)
	{
		master = make_shared<SecureRandom>(make_shared<CryptoApiRandomGenerator>());
	}
	DoubleScale = pow(2.0, 64.0);
}

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
//#include "../Interfaces/ISecureRandom.h"

using namespace std;


class BigInteger
{
private:
	// Resource string
	static const char * DivisionByZero;
	static const char * ModulusPositive;
	static const char * NotRelativelyPrime;
	static const char * NegativeValue;
	static const char * NegativeExponent;
	static const char * ResultTooLarge;
	static const char * NegativeBitPosition;
	static const char * InvalidBitAddress;
	static const char * ZeroLengthBigInteger;
	static const char * InvalidSign;
	static const char * NegativeSizeInBits;
	static const char * InvalidBitLength;
	static const char * InvalidBase;
	static const char * BadCharacterRadix8;
	static const char * BadCharacterRadix2;
	static const char * UnSupportedBase;

private:
	static const int64_t IMASK = 0xFFFFFFFF;
	static const uint64_t UIMASK = 0xFFFFFFFF;
	static vector<uint8_t> BitLengthTable;

	/// <summary>
	/// These are the threshold bit-lengths (of an exponent) where we
	/// increase the window size. <br />They are calculated according to the
	/// expected savings in multiplications. <br />Some squares will also be
	/// saved on average, but we offset these against the extra storage
	/// costs. <br />
	/// </summary>

	static vector<int32_t> ExpWindowThresholds;

	// TODO Parse radix-2 64 bits at a time and radix-8 63 bits at a time
	static const int32_t chunk2 = 1;
	static const int32_t chunk8 = 1;
	static const int32_t chunk10 = 19;
	static const int32_t chunk16 = 16;

	static const int32_t BitsPerByte = 8;
	static const int32_t BitsPerInt = 32;
	static const int32_t BytesPerInt = 4;

private:
	// array of ints with [0] being the most significant
	vector<int32_t> magnitude;
	int32_t sign; // -1 means -ve; +1 means +ve; 0 means 0;
	int32_t nBits; // cache BitCount() value
	int32_t nBitLength; // cache BitLength() value
						// -m^(-1) mod b, b = 2^32 (see Montgomery mult.), 0 when uninitialised
	int32_t mQuote;
	bool IsInitialized = false;

private:
	static BigInteger Zero, One, Two, Three, Ten;
	// Each list has a product < 2^31
	static vector<vector<int32_t>> primeLists;
	static vector<int32_t> primeProducts, ZeroMagnitude;
	static vector<uint8_t> ZeroEncoding;
	static vector<BigInteger> SMALL_CONSTANTS;
	static BigInteger radix2, radix2E, radix8, radix8E, radix10, radix10E, radix16, radix16E;
	static ISecureRandom RandomSource;


private:
	BigInteger AddToMagnitude(const vector<int32_t> &magToAdd) const;
	inline bool QuickPow2Check() const;

	/// <summary>
	/// return z = x / y - done in place (z value preserved, x contains the *
	/// remainder)
	/// </summary>
	vector<int32_t> Divide(vector<int32_t> &x, vector<int32_t> &y) const;
	bool IsEqualMagnitude(const BigInteger &x) const;

	bool CheckProbablePrime(const int32_t certainty, IRandom random, const bool randomlySelected);

	BigInteger ModInversePow2(const BigInteger &value) const;

	/// <summary>
	/// Calculate mQuote = -m^(-1) mod b with b = 2^32 (32 = word size)
	/// </summary>
	inline int32_t GetMQuote();

	inline int32_t Remainder(const int32_t m) const;

	/// <summary>
	/// return x = x mod y - done in place (y value preserved)
	/// </summary>
	vector<int32_t> Remainder(vector<int32_t> &x, const vector<int32_t> &y) const;

	inline vector<int32_t> LastNBits(const int32_t n) const;

	inline BigInteger DivideWords(const int32_t w) const;

	inline BigInteger RemainderWords(const int32_t w) const;

	vector<uint8_t> ToByteArray(const bool _unsigned) const;

	inline BigInteger FlipExistingBit(const int32_t n) const;

	int32_t GetLowestSetBitMaskFirst(const int32_t firstWordMask) const;

	void ParseString(const string &str, const int32_t radix);
	void ParseBytes(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length);
	void ParseBytesWithSign(const int32_t sign, const vector<uint8_t> &bytes, const int32_t offset, const int32_t length);

	inline static int32_t GetByteLength(const int32_t nBits);

	static vector<int32_t> MakeMagnitude(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length);

	/// <summary>
	/// a = a + b - b preserved.
	/// </summary>
	inline static vector<int32_t> AddMagnitudes(vector<int32_t> &a, const vector<int32_t> &b);

	static int32_t CalcBitLength(const int32_t sign, const int32_t indx, const vector<int32_t> &mag);

	/// <summary>
	/// unsigned comparison on two arrays - note the arrays may start with
	/// leading zeros.
	/// </summary>
	static int32_t CompareTo(const int32_t xIndx, const vector<int32_t> &x, const int32_t yIndx, const vector<int32_t> &y);

	static int32_t CompareNoLeadingZeroes(const int32_t xIndx, const vector<int32_t> &x, const int32_t yIndx, const vector<int32_t> &y);

	inline static int32_t ModInverse32(const int32_t d);

	inline static int64_t ModInverse64(const int64_t d);

	/// <summary>
	/// Calculate the numbers u1, u2, and u3 such that: <br />u1 * a + u2 * b
	/// = u3 <br />where u3 is the greatest common divider of a and b. <br />
	/// a and b using the extended Euclid algorithm (refer p. 323 of The Art
	/// of Computer Programming vol 2, 2nd ed). <br />This also seems to have
	/// the side effect of calculating some form of multiplicative inverse.
	/// </summary>
	/// <param name="a">
	/// First number to calculate gcd for
	/// </param>
	/// <param name="b">
	/// Second number to calculate gcd for
	/// </param>
	/// <param name="u1Out">
	/// the return object for the u1 value
	/// </param>
	/// <returns>
	/// The greatest common divisor of a and b
	/// </returns>
	inline static BigInteger ExtEuclid(const BigInteger &a, const BigInteger &b, BigInteger &u1Out);

	static BigInteger ModPowBarrett(const BigInteger &b, const BigInteger &e, const BigInteger &m);

	static BigInteger ReduceBarrett(BigInteger &x, BigInteger &m, const BigInteger &mr, const BigInteger &yu);

	static BigInteger ModPowMonty(BigInteger &b, const BigInteger &e, const BigInteger &m, const bool convert);

	static vector<int32_t> GetWindowList(const vector<int32_t> &mag, const int32_t extraBits);

	inline static int32_t CreateWindowEntry(int32_t mult, int32_t zeroes);

	/// <returns>
	/// w with w = x * x - w is assumed to have enough space.
	/// </returns>
	static vector<int32_t> Square(vector<int32_t> &w, const vector<int32_t> &x);

	/// <returns>
	/// x with x = y * z - x is assumed to have enough space.
	/// </returns>
	static vector<int32_t> Multiply(vector<int32_t> &x, const vector<int32_t> &y, const vector<int32_t> &z);

	// mDash = -m^(-1) mod b
	static void MontgomeryReduce(vector<int32_t> &x, const vector<int32_t> &m, const uint32_t mDash);

	// mDash = -m^(-1) mod b

	/// <summary>
	/// Montgomery multiplication: a = x * y * R^(-1) mod m <br />Based
	/// algorithm 14.36 of Handbook of Applied Cryptography. <br />&lt;li&gt;
	/// m, x, y should have length n &lt;/li&gt; <br />&lt;li&gt; a should
	/// have length (n + 1) &lt;/li&gt; <br />&lt;li&gt; b = 2^32, R = b^n
	/// &lt;/li&gt; <br />&lt;br/&gt; <br />The result is put in x <br />
	/// &lt;br/&gt; <br />NOTE: the indices of x, y, m, a different in HAC
	/// and in Java <br />
	/// </summary>
	static void MultiplyMonty(vector<int32_t> &a, vector<int32_t> &x, const vector<int32_t> &y,
		const vector<int32_t> &m, const uint32_t mDash, const bool smallMontyModulus);

	// mDash = -m^(-1) mod b
	static void SquareMonty(vector<int32_t> &a, vector<int32_t> &x, const vector<int32_t> &m,
		const uint32_t mDash, const bool smallMontyModulus);

	inline static uint32_t MultiplyMontyNIsOne(const uint32_t x, const uint32_t y, const uint32_t m, const uint32_t mDash);

	/// <summary>
	/// do a left shift - this returns a new array.
	/// </summary>
	static vector<int32_t> ShiftLeft(const vector<int32_t> &mag, const int32_t n);

	/// <summary>
	/// do a right shift - this does it in place.
	/// </summary>
	static void ShiftRightInPlace(const int32_t start, vector<int32_t> &mag, const int32_t n);

	/// <summary>
	/// do a right shift by one - this does it in place.
	/// </summary>
	static void ShiftRightOneInPlace(const int32_t start, vector<int32_t> &mag);

	/// <summary>
	/// returns x = x - y - we assume x is &gt;= y
	/// </summary>
	static vector<int32_t> Subtract(const int32_t xStart, vector<int32_t> &x, const int32_t yStart, const vector<int32_t> &y);

	inline static vector<int32_t> doSubBigLil(const vector<int32_t> &bigMag, const vector<int32_t> &lilMag);

	inline static void AppendZeroExtendedString(ostringstream &sl, const string s, const int32_t minLength);

	static void ToString(ostringstream &sl, int32_t radix, vector<BigInteger> &moduli, int32_t scale, BigInteger &pos);

	static BigInteger CreateUValueOf(const uint64_t value);
	static BigInteger CreateValueOf(const int64_t value);

	static string IntToBin(const int32_t input);

	static string IntToOctal(const int32_t input);

private:
	BigInteger(const int32_t signum, const vector<int32_t> &mag, const bool checkMag);

public:
	BigInteger();
	BigInteger(const string &value);
	BigInteger(const string &value, const int32_t radix);
	BigInteger(const vector<uint8_t> &bytes);
	BigInteger(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length);
	BigInteger(const int32_t sign, const vector<uint8_t> &bytes);
	BigInteger(const int32_t sign, const vector<uint8_t> &bytes, const int32_t offset, const int32_t length);
	BigInteger(const int32_t sizeInBits, IRandom random);
	BigInteger(const int32_t BitLength, const int32_t certainty, IRandom random);

	BigInteger Abs() const;
	BigInteger Add(const BigInteger &value) const;
	BigInteger Subtract(const BigInteger &n) const;
	BigInteger And(const BigInteger &value) const;
	BigInteger Not() const;
	BigInteger AndNot(const BigInteger &val) const;
	BigInteger Or(const BigInteger &value) const;
	BigInteger Xor(const BigInteger &value) const;

	int32_t CompareTo(const BigInteger &value) const;
	BigInteger Divide(const BigInteger &val) const;
	vector<BigInteger> DivideAndRemainder(const BigInteger &val) const;
	BigInteger Gcd(const BigInteger &value) const;
	BigInteger Inc() const;

	bool RabinMillerTest(const int32_t certainty, IRandom random) const;

	bool RabinMillerTest(const int32_t certainty, IRandom random, const bool randomlySelected) const;

	/// <summary>
	/// return whether or not a BigInteger is probably prime with a
	/// probability of 1 - (1/2)**certainty. <br />&lt;p&gt;From Knuth Vol 2,
	/// pg 395.&lt;/p&gt;
	/// </summary>
	bool IsProbablePrime(const int32_t certainty) const;
	bool IsProbablePrime(const int32_t certainty, const bool randomlySelected) const;

	BigInteger Max(const BigInteger &value) const;
	BigInteger Min(const BigInteger &value) const;
	BigInteger Mod(const BigInteger &m) const;
	BigInteger ModInverse(const BigInteger &m) const;
	BigInteger ModPow(BigInteger &e, const BigInteger &m) const;

	BigInteger Multiply(const BigInteger &val) const;
	BigInteger Square() const;
	BigInteger Negate() const;

	BigInteger NextProbablePrime() const;

	BigInteger Pow(const int32_t exp) const;

	BigInteger Remainder(const BigInteger &n) const;

	BigInteger ShiftLeft(const int32_t n) const;
	BigInteger ShiftRight(const int32_t n) const;

	vector<uint8_t> ToByteArray() const;
	vector<uint8_t> ToByteArrayUnsigned() const;

	bool TestBit(const int32_t n) const;
	BigInteger SetBit(const int32_t n) const;
	BigInteger ClearBit(const int32_t n) const;
	BigInteger FlipBit(const int32_t n) const;

	int32_t GetLowestSetBit() const;

	string ToString() const;
	string ToString(const int32_t radix) const;

	inline bool Equals(const BigInteger &other) const;
	int32_t GetHashCode() const;
	static int32_t BitCnt(const int32_t i);

	static int32_t BitLen(const int32_t w);
	static BigInteger ProbablePrime(const int32_t BitLength, IRandom random);

	static BigInteger ValueOf(const int64_t value);

	static BigInteger Arbitrary(const int32_t sizeInBits);

	static void Boot();

public:
	inline int32_t GetBitLength();
	inline int32_t GetBitCount();
	inline int32_t GetInt32Value() const;
	inline int64_t GetInt64Value() const;

	inline bool GetIsInitialized() const
	{
		return IsInitialized;
	}

	inline int32_t GetSignValue() const
	{
		return sign;
	}

	inline static BigInteger GetZero()
	{
		return Zero;
	}

	inline static BigInteger GetOne()
	{
		return One;
	}

	inline static BigInteger GetTwo()
	{
		return Two;
	}

	inline static BigInteger GetThree()
	{
		return Three;
	}

	inline static BigInteger GetTen()
	{
		return Ten;
	}

	inline static vector<vector<int32_t>> GetprimeLists()
	{
		return primeLists;
	}

	inline static vector<int32_t> GetprimeProducts()
	{
		return primeProducts;
	}

	inline static IRandom GetRandomSource()
	{
		return RandomSource;
	}

};


// ///////////////////////////////////////////////////////////////// //
// *C++ 11 BigInteger Library                                 
// *Copyright(c) 2018  Mbadiwe Nnaemeka Ronald                 
// *Github Repository <https://github.com/ron4fun>             

// *Distributed under the MIT software license, see the accompanying file LICENSE 
// *or visit http ://www.opensource.org/licenses/mit-license.php.           

// *Acknowledgements:                                  
// ** //
// *Thanks to Ugochukwu Mmaduekwe (https://github.com/Xor-el) for his creative        
// *development of this library in Pascal/Delphi                         

// ////////////////////////////////////////////////////// ///////////////
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

const char * BigInteger::DivisionByZero = "Division by Zero Error";
const char * BigInteger::ModulusPositive = "Modulus must be Positive";
const char * BigInteger::NotRelativelyPrime = "Numbers not Relatively Prime";
const char * BigInteger::NegativeValue = "Cannot be Called on Value < 0";
const char * BigInteger::NegativeExponent = "Negative Exponent";
const char * BigInteger::ResultTooLarge = "Result too Large";
const char * BigInteger::NegativeBitPosition = "Bit Position must not be Negative";
const char * BigInteger::InvalidBitAddress = "Bit Address less than Zero";
const char * BigInteger::ZeroLengthBigInteger = "Zero length BigInteger";
const char * BigInteger::InvalidSign = "Invalid Sign Value";
const char * BigInteger::NegativeSizeInBits = "sizeInBits must be non-negative";
const char * BigInteger::InvalidBitLength = "bitLength < 2";
const char * BigInteger::InvalidBase = "Only bases 2, 8, 10, or 16 allowed";
const char * BigInteger::BadCharacterRadix8 = "Bad Character in radix 8 string: %s";
const char * BigInteger::BadCharacterRadix2 = "Bad Character in radix 2 string: %s";
const char * BigInteger::UnSupportedBase = "Only bases 2, 8, 10, 16 are allowed";

BigInteger BigInteger::Zero = BigInteger(), BigInteger::One = BigInteger(), BigInteger::Two = BigInteger(),
BigInteger::Three = BigInteger(), BigInteger::Ten = BigInteger();

vector<vector<int32_t>> BigInteger::primeLists = vector<vector<int32_t>>();

vector<int32_t> BigInteger::primeProducts = vector<int32_t>();

vector<int32_t> BigInteger::ZeroMagnitude = vector<int32_t>();

vector<uint8_t> BigInteger::ZeroEncoding = vector<uint8_t>();

vector<BigInteger> BigInteger::SMALL_CONSTANTS = vector<BigInteger>();

BigInteger BigInteger::radix2 = BigInteger(), BigInteger::radix2E = BigInteger(), BigInteger::radix8 = BigInteger(),
BigInteger::radix8E = BigInteger(), BigInteger::radix10 = BigInteger(), BigInteger::radix10E = BigInteger(),
BigInteger::radix16 = BigInteger(), BigInteger::radix16E = BigInteger();

ISecureRandom BigInteger::RandomSource = nullptr;

vector<uint8_t> BigInteger::BitLengthTable = vector<uint8_t>({ 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4,
	4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
	6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 });

vector<int32_t> BigInteger::ExpWindowThresholds = vector<int32_t>({ 7, 25, 81, 241, 673, 1793, 4609, int32_t(0x7FFFFFFF) });

BigInteger::BigInteger()
{
	// Boot
}

BigInteger::BigInteger(const string &value)
{
	ParseString(value, 10);
}

BigInteger::BigInteger(const string &value, const int32_t radix)
{
	ParseString(value, radix);
}

BigInteger::BigInteger(const vector<uint8_t> &bytes)
{
	ParseBytes(bytes, 0, bytes.size());
}

BigInteger::BigInteger(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length)
{
	ParseBytes(bytes, offset, length);
}

BigInteger::BigInteger(const int32_t sign, const vector<uint8_t> &bytes)
{
	ParseBytesWithSign(sign, bytes, 0, bytes.size());
}

BigInteger::BigInteger(const int32_t sign, const vector<uint8_t> &bytes, const int32_t offset, const int32_t length)
{
	ParseBytesWithSign(sign, bytes, offset, length);
}

BigInteger::BigInteger(const int32_t sizeInBits, IRandom random)
{
	int32_t nBytes, xBits;
	vector<uint8_t> b;

	if (sizeInBits < 0)
		throw invalid_argument(NegativeSizeInBits);

	sign = -1;
	nBits = -1;
	nBitLength = -1;
	mQuote = 0;
	IsInitialized = true;

	if (sizeInBits == 0)
	{
		sign = 0;
		magnitude = ZeroMagnitude;
		return;
	}

	nBytes = GetByteLength(sizeInBits);
	b = vector<uint8_t>(nBytes);
	random.get()->NextBytes(b);

	// strip off any excess bits in the MSB
	xBits = (BitsPerByte * nBytes) - sizeInBits;
	b[0] = b[0] & uint8_t(uint32_t(255) >> xBits);

	magnitude = MakeMagnitude(b, 0, b.size());

	if (magnitude.size() < 1) sign = 0;
	else
		sign = 1;
}

BigInteger::BigInteger(const int32_t BitLength, const int32_t certainty, IRandom random)
{
	int32_t nBytes, xBits, j;
	uint8_t mask, lead;
	vector<uint8_t> b;

	if (BitLength < 2)
		throw invalid_argument(InvalidBitLength);

	sign = 1;
	nBits = -1;
	nBitLength = BitLength;
	mQuote = 0;
	IsInitialized = true;

	if (BitLength == 2)
	{
		if (random.get()->Next(2) == 0) magnitude = Two.magnitude;
		else
			magnitude = Three.magnitude;
	}

	nBytes = GetByteLength(BitLength);
	b = vector<uint8_t>(nBytes);

	xBits = (BitsPerByte * nBytes) - BitLength;
	mask = uint8_t(uint32_t(255) >> xBits);
	lead = uint8_t(1 << (7 - xBits));

	while (true)
	{
		random.get()->NextBytes(b);

		// strip off any excess bits in the MSB
		b[0] = b[0] & mask;

		// ensure the leading bit is 1 (to meet the strength requirement)
		b[0] = b[0] | lead;

		// ensure the trailing bit is 1 (i.e. must be odd)
		b[nBytes - 1] = b[nBytes - 1] | 1;

		magnitude = MakeMagnitude(b, 0, b.size());
		nBits = -1;
		mQuote = 0;

		if (certainty < 1) break;

		if (CheckProbablePrime(certainty, random, true)) break;

		j = 1;
		while (j < (magnitude.size() - 1))
		{
			magnitude[j] = magnitude[j] ^ random.get()->Next();

			if (CheckProbablePrime(certainty, random, true)) return;

			j++;
		}

	}

}

BigInteger::BigInteger(const int32_t signum, const vector<int32_t> &mag, const bool checkMag)
{
	int32_t i;

	sign = -1;
	nBits = -1;
	nBitLength = -1;
	mQuote = 0;
	IsInitialized = true;

	if (checkMag)
	{
		i = 0;
		while ((i < mag.size()) && (mag[i] == 0)) i++;

		if (i == mag.size())
		{
			sign = 0;
			magnitude = ZeroMagnitude;
		}
		else
		{
			sign = signum;

			if (i == 0) magnitude = mag;
			else
			{
				// strip leading 0 words
				magnitude = vector<int32_t>(mag.size() - i);
				if (!magnitude.empty())
					memmove(&magnitude[0], &mag[i], magnitude.size() * sizeof(int32_t));
			}

		}

	}
	else
	{
		sign = signum;
		magnitude = mag;
	}

}


BigInteger BigInteger::Abs() const
{
	if (sign >= 0) return *this;

	return Negate();
}

BigInteger BigInteger::Add(const BigInteger &value) const
{
	if (sign == 0) return value;

	if (sign != value.sign)
	{
		if (value.sign == 0) return *this;

		if (value.sign < 0) return Subtract(value.Negate());

		return value.Subtract(Negate());
	}

	return AddToMagnitude(value.magnitude);
}

BigInteger BigInteger::Subtract(const BigInteger &n) const
{
	int32_t compare;
	BigInteger bigun, lilun;

	if (n.sign == 0) return *this;

	if (sign == 0) return n.Negate();

	if (sign != n.sign) return Add(n.Negate());

	compare = CompareNoLeadingZeroes(0, magnitude, 0, n.magnitude);
	if (compare == 0) return Zero;

	if (compare < 0)
	{
		bigun = n;
		lilun = *this;
	}
	else
	{
		bigun = *this;
		lilun = n;
	}

	return BigInteger(sign * compare, doSubBigLil(bigun.magnitude, lilun.magnitude), true);
}

BigInteger BigInteger::And(const BigInteger &value) const
{
	vector<int32_t> aMag, bMag, resultMag;
	bool resultNeg;
	int32_t resultLength, aStart, bStart, i, aWord, bWord;

	if ((sign == 0) || (value.sign == 0)) return Zero;

	if (sign > 0) aMag = magnitude;
	else
		aMag = Add(One).magnitude;

	if (value.sign > 0) bMag = value.magnitude;
	else
		bMag = value.Add(One).magnitude;

	resultNeg = (sign < 0) && (value.sign < 0);
	resultLength = std::max(aMag.size(), bMag.size());

	resultMag = vector<int32_t>(resultLength);

	aStart = resultMag.size() - aMag.size();
	bStart = resultMag.size() - bMag.size();

	for (i = 0; i < resultMag.size(); i++)
	{
		if (i >= aStart) aWord = aMag[i - aStart];
		else
			aWord = 0;

		if (i >= bStart) bWord = bMag[i - bStart];
		else
			bWord = 0;

		if (sign < 0) aWord = ~aWord;

		if (value.sign < 0) bWord = ~bWord;

		resultMag[i] = aWord & bWord;

		if (resultNeg) resultMag[i] = ~resultMag[i];
	}

	BigInteger result = BigInteger(1, resultMag, true);

	// TODO Optimise this case
	if (resultNeg) result = result.Not();

	return result;
}

BigInteger BigInteger::Not() const
{
	return Inc().Negate();
}

BigInteger BigInteger::AndNot(const BigInteger &val) const
{
	return And(val.Not());
}

BigInteger BigInteger::Or(const BigInteger &value) const
{
	vector<int32_t> aMag, bMag, resultMag;
	bool resultNeg;
	int32_t resultLength, aStart, bStart, i, aWord, bWord;

	if (sign == 0) return value;

	if (value.sign == 0) return *this;

	if (sign > 0) aMag = magnitude;
	else
		aMag = Add(One).magnitude;

	if (value.sign > 0) bMag = value.magnitude;
	else
		bMag = value.Add(One).magnitude;

	resultNeg = (sign < 0) || (value.sign < 0);
	resultLength = std::max(aMag.size(), bMag.size());

	resultMag = vector<int32_t>(resultLength);

	aStart = resultMag.size() - aMag.size();
	bStart = resultMag.size() - bMag.size();

	for (i = 0; i < resultMag.size(); i++)
	{
		if (i >= aStart) aWord = aMag[i - aStart];
		else
			aWord = 0;

		if (i >= bStart) bWord = bMag[i - bStart];
		else
			bWord = 0;

		if (sign < 0) aWord = ~aWord;

		if (value.sign < 0) bWord = ~bWord;

		resultMag[i] = aWord | bWord;

		if (resultNeg) resultMag[i] = ~resultMag[i];
	}

	BigInteger result = BigInteger(1, resultMag, true);

	// TODO Optimise this case
	if (resultNeg) result = result.Not();

	return result;
}

BigInteger BigInteger::Xor(const BigInteger &value) const
{
	vector<int32_t> aMag, bMag, resultMag;
	bool resultNeg;
	int32_t resultLength, aStart, bStart, i, aWord, bWord;

	if (sign == 0) return value;

	if (value.sign == 0) return *this;

	if (sign > 0) aMag = magnitude;
	else
		aMag = Add(One).magnitude;

	if (value.sign > 0) bMag = value.magnitude;
	else
		bMag = value.Add(One).magnitude;

	// TODO Can just replace with sign != value.sign?
	resultNeg = ((sign < 0) && (value.sign >= 0)) || ((sign >= 0) && (value.sign < 0));
	resultLength = std::max(aMag.size(), bMag.size());

	resultMag = vector<int32_t>(resultLength);

	aStart = resultMag.size() - aMag.size();
	bStart = resultMag.size() - bMag.size();

	for (i = 0; i < resultMag.size(); i++)
	{
		if (i >= aStart) aWord = aMag[i - aStart];
		else
			aWord = 0;

		if (i >= bStart) bWord = bMag[i - bStart];
		else
			bWord = 0;

		if (sign < 0) aWord = ~aWord;

		if (value.sign < 0) bWord = ~bWord;

		resultMag[i] = aWord ^ bWord;

		if (resultNeg) resultMag[i] = ~resultMag[i];
	}

	BigInteger result = BigInteger(1, resultMag, true);

	// TODO Optimise this case
	if (resultNeg) result = result.Not();

	return result;
}

int32_t BigInteger::CompareTo(const BigInteger &value) const
{
	if (sign < value.sign) return -1;

	if (sign > value.sign) return 1;

	if (sign == 0) return 0;

	return sign * CompareNoLeadingZeroes(0, magnitude, 0, value.magnitude);
}

BigInteger BigInteger::Divide(const BigInteger &val) const
{
	BigInteger tempRes;
	vector<int32_t> mag, valMag;

	if (val.sign == 0)
		throw invalid_argument(DivisionByZero);

	if (sign == 0) return Zero;

	if (val.QuickPow2Check()) // val is power of two
	{
		tempRes = Abs().ShiftRight(val.Abs().GetBitLength() - 1);
		if (val.sign == sign) return tempRes;

		return tempRes.Negate();
	}

	mag = magnitude;
	valMag = val.magnitude;
	return BigInteger(sign * val.sign, Divide(mag, valMag), true);
}

vector<BigInteger> BigInteger::DivideAndRemainder(const BigInteger &val) const
{
	vector<BigInteger> biggies;
	int32_t e;
	BigInteger Quotient;
	vector<int32_t> Remainder, quotient_array, valMag;

	if (val.sign == 0)
		throw invalid_argument(DivisionByZero);

	biggies = vector<BigInteger>(2);

	if (sign == 0)
	{
		biggies[0] = Zero;
		biggies[1] = Zero;
	}
	else if (val.QuickPow2Check()) // val is power of two
	{
		e = val.Abs().GetBitLength() - 1;
		Quotient = Abs().ShiftRight(e);
		Remainder = LastNBits(e);

		if (val.sign == sign) biggies[0] = Quotient;
		else
			biggies[0] = Quotient.Negate();

		biggies[1] = BigInteger(sign, Remainder, true);
	}
	else
	{
		Remainder = magnitude;
		valMag = val.magnitude;
		quotient_array = Divide(Remainder, valMag);

		biggies[0] = BigInteger(sign * val.sign, quotient_array, true);
		biggies[1] = BigInteger(sign, Remainder, true);
	}

	return biggies;
}

BigInteger BigInteger::Gcd(const BigInteger &value) const
{
	BigInteger r, u, v;

	if (value.sign == 0) return Abs();

	if (sign == 0) return value.Abs();

	u = *this;
	v = value;

	while (v.sign != 0)
	{
		r = u.Mod(v);
		u = v;
		v = r;
	}

	return u;
}

BigInteger BigInteger::Inc() const
{
	if (sign == 0) return One;

	if (sign < 0)
		return BigInteger(-1, doSubBigLil(magnitude, One.magnitude), true);

	return AddToMagnitude(One.magnitude);
}

bool BigInteger::RabinMillerTest(const int32_t certainty, IRandom random) const
{
	return RabinMillerTest(certainty, random, false);
}

bool BigInteger::RabinMillerTest(const int32_t certainty, IRandom random, const bool randomlySelected) const
{
	int32_t bits, iterations, itersFor100Cert, s, j, shiftval;
	BigInteger r, montRadix, minusMontRadix, a, y, temp = *this;

	bits = temp.GetBitLength();

	iterations = ((certainty - 1) / 2) + 1;
	if (randomlySelected)
	{
		if (bits >= 1024) itersFor100Cert = 4;
		else if (bits >= 512) itersFor100Cert = 8;
		else if (bits >= 256) itersFor100Cert = 16;
		else
			itersFor100Cert = 50;

		if (certainty < 100)
			iterations = std::min(itersFor100Cert, iterations);
		else
		{
			iterations = iterations - 50;
			iterations = iterations + itersFor100Cert;
		}

	}

	// let n = 1 + d . 2^s
	shiftval = int32_t(-1) << 1; // -2
	s = temp.GetLowestSetBitMaskFirst(shiftval);
	r = temp.ShiftRight(s);

	// NOTE: Avoid conversion to/from Montgomery form and check for R/-R as result instead

	montRadix = One.ShiftLeft(32 * temp.magnitude.size()).Remainder(temp);
	minusMontRadix = temp.Subtract(montRadix);

	do
	{
		do
		{
			a = BigInteger(temp.GetBitLength(), random);

		} while (((a.sign == 0) || (a.CompareTo(temp) >= 0) ||
			(a.IsEqualMagnitude(montRadix)) || (a.IsEqualMagnitude(minusMontRadix))));

		y = ModPowMonty(a, r, temp, false);

		if (!y.Equals(montRadix))
		{
			j = 0;
			while (!y.Equals(minusMontRadix))
			{
				j++;
				if (j == s) return false;

				y = ModPowMonty(y, Two, temp, false);

				if (y.Equals(montRadix)) return false;
			}
		}

		iterations--;
	} while (iterations > 0);

	return true;
}

bool BigInteger::IsProbablePrime(const int32_t certainty) const
{
	return IsProbablePrime(certainty, false);
}

bool BigInteger::IsProbablePrime(const int32_t certainty, const bool randomlySelected) const
{
	BigInteger n;

	if (certainty <= 0) return true;

	n = Abs();

	if (!n.TestBit(0)) return n.Equals(Two);

	if (n.Equals(One)) return false;

	return n.CheckProbablePrime(certainty, RandomSource, randomlySelected);
}

BigInteger BigInteger::Max(const BigInteger &value) const
{
	if (CompareTo(value) > 0) return *this;

	return value;
}

BigInteger BigInteger::Min(const BigInteger &value) const
{
	if (CompareTo(value) < 0) return *this;

	return value;
}

BigInteger BigInteger::Mod(const BigInteger &m) const
{
	BigInteger biggie;

	if (m.sign < 1)
		throw invalid_argument(ModulusPositive);

	biggie = Remainder(m);

	if (biggie.sign >= 0) return biggie;

	return biggie.Add(m);
}

BigInteger BigInteger::ModInverse(const BigInteger &m) const
{
	BigInteger d, x, Gcd;

	if (m.sign < 1)
		throw invalid_argument(ModulusPositive);

	if (m.QuickPow2Check()) return ModInversePow2(m);

	d = Remainder(m);
	Gcd = ExtEuclid(d, m, x);

	if (!Gcd.Equals(One))
		throw invalid_argument(NotRelativelyPrime);

	if (x.sign < 0) x = x.Add(m);

	return x;
}

BigInteger BigInteger::ModPow(BigInteger &e, const BigInteger &m) const
{
	bool negExp;

	if (m.sign < 1)
		throw invalid_argument(ModulusPositive);

	if (m.Equals(One)) return Zero;

	if (e.sign == 0) return One;

	if (sign == 0) return Zero;

	negExp = e.sign < 0;
	if (negExp) e = e.Negate();

	BigInteger result = Mod(m);

	if (!e.Equals(One))
	{
		if ((m.magnitude[m.magnitude.size() - 1] & 1) == 0)
			result = ModPowBarrett(result, e, m);
		else
			result = ModPowMonty(result, e, m, true);
	}

	if (negExp) result = result.ModInverse(m);

	return result;
}

BigInteger BigInteger::Multiply(const BigInteger &val) const
{
	int32_t resLength, resSign;
	vector<int32_t> res;

	if (val.Equals(*this)) return Square();

	if ((sign & val.sign) == 0) return Zero;

	BigInteger result;
	if (val.QuickPow2Check()) // val is power of two
	{
		result = ShiftLeft(val.Abs().GetBitLength() - 1);
		if (val.sign > 0) return result;

		return result.Negate();
	}

	if (QuickPow2Check()) // this is power of two
	{
		result = val.ShiftLeft(Abs().GetBitLength() - 1);
		if (sign > 0) return result;

		return result.Negate();
	}

	resLength = magnitude.size() + val.magnitude.size();
	res = vector<int32_t>(resLength);

	Multiply(res, magnitude, val.magnitude);

	resSign = sign ^ val.sign ^ 1;

	return BigInteger(resSign, res, true);
}

BigInteger BigInteger::Square() const
{
	int32_t resLength;
	vector<int32_t> res;

	if (sign == 0) return Zero;

	if (QuickPow2Check()) return ShiftLeft(Abs().GetBitLength() - 1);

	resLength = magnitude.size() << 1;

	if (uint32_t(magnitude[0]) >> 16 == 0) resLength--;

	res = vector<int32_t>(resLength);
	Square(res, magnitude);

	return BigInteger(1, res, false);
}

BigInteger BigInteger::Negate() const
{
	if (sign == 0) return *this;

	return BigInteger(-sign, magnitude, false);
}

BigInteger BigInteger::NextProbablePrime() const
{
	BigInteger n;

	if (sign < 0)
		throw invalid_argument(NegativeValue);

	if (CompareTo(Two) < 0) return Two;

	n = Inc().SetBit(0);

	while (!n.CheckProbablePrime(100, RandomSource, false))
		n = n.Add(Two);

	return n;
}

BigInteger BigInteger::Pow(const int32_t _exp) const
{
	int64_t powOf2;
	BigInteger y, z;
	int32_t exp = _exp;

	if (exp <= 0)
	{
		if (exp < 0)
			throw invalid_argument(NegativeExponent);

		return One;
	}

	if (sign == 0) return *this;

	z = *this;
	if (QuickPow2Check())
	{
		powOf2 = int64_t(exp) * (z.GetBitLength() - 1);
		if (powOf2 > INT32_MAX)
			throw invalid_argument(ResultTooLarge);

		return One.ShiftLeft(int32_t(powOf2));
	}

	y = One;
	while (true)
	{
		if ((exp & 0x1) == 1) y = y.Multiply(z);

		exp = exp >> 1;
		if (exp == 0) break;

		z = z.Multiply(z);
	}

	return y;
}

BigInteger BigInteger::Remainder(const BigInteger &n) const
{
	int32_t val, rem;
	vector<int32_t> tempRes;

	if (n.sign == 0)
		throw invalid_argument(DivisionByZero);

	if (sign == 0) return Zero;

	// For small values, use fast remainder method
	if (n.magnitude.size() == 1)
	{
		val = n.magnitude[0];

		if (val > 0)
		{
			if (val == 1) return Zero;

			// TODO Make this func work on uint, and handle val == 1?
			rem = Remainder(val);

			if (rem == 0) return Zero;

			return BigInteger(sign, vector<int32_t>({ rem }), false);
		}
	}

	if (CompareNoLeadingZeroes(0, magnitude, 0, n.magnitude) < 0)
		return *this;

	if (n.QuickPow2Check()) // n is power of two
							// TODO Move before small values branch above?
		tempRes = LastNBits(n.Abs().GetBitLength() - 1);
	else
	{
		tempRes = magnitude;
		tempRes = Remainder(tempRes, n.magnitude);
	}

	return BigInteger(sign, tempRes, true);
}

int32_t BigInteger::Remainder(const int32_t m) const
{
	int64_t acc, posVal;
	int32_t pos;

	acc = 0;
	for (pos = 0; pos < magnitude.size(); pos++)
	{
		posVal = uint32_t(magnitude[pos]);
		acc = ((uint64_t(acc) << 32) | posVal) % m;
	}

	return int32_t(acc);
}

BigInteger BigInteger::ShiftLeft(const int32_t n) const
{
	if ((sign == 0) || (magnitude.size() == 0))
		return Zero;

	if (n == 0) return *this;

	if (n < 0) return ShiftRight(-n);

	BigInteger result = BigInteger(sign, ShiftLeft(magnitude, n), true);
	BigInteger temp = *this;

	if (temp.GetBitCount() != -1)
	{
		if (sign > 0) result.nBits = temp.nBits;
		else
			result.nBits = temp.nBits + n;
	}

	if (temp.nBitLength != -1) result.nBitLength = temp.nBitLength + n;

	return result;
}

BigInteger BigInteger::ShiftRight(const int32_t n) const
{
	int32_t resultLength, numInts, numBits, numBits2, magPos, i;
	vector<int32_t>	res;

	BigInteger temp = *this;
	if (n == 0) return temp;

	if (n < 0) return ShiftLeft(-n);

	if (n >= temp.GetBitLength())
	{
		if (sign < 0) return One.Negate();

		return Zero;
	}

	resultLength = (temp.GetBitLength() - n + 31) >> 5;
	res = vector<int32_t>(resultLength);

	numInts = int32_t(uint32_t(n) >> 5);
	numBits = n & 31;

	if (numBits == 0)
		memmove(&res[0], &magnitude[0], res.size() * sizeof(int32_t));
	else
	{
		numBits2 = 32 - numBits;

		magPos = magnitude.size() - 1 - numInts;

		i = resultLength - 1;
		while (i >= 0)
		{
			res[i] = int32_t(uint32_t(magnitude[magPos]) >> numBits);
			magPos--;

			if (magPos >= 0)
				res[i] = res[i] | (uint32_t(magnitude[magPos]) << numBits2);

			i--;
		}
	}

	return BigInteger(sign, res, false);
}

vector<uint8_t> BigInteger::ToByteArray() const
{
	return ToByteArray(false);
}

vector<uint8_t> BigInteger::ToByteArrayUnsigned() const
{
	return ToByteArray(true);
}

bool BigInteger::TestBit(const int32_t n) const
{
	int32_t wordNum, word;

	if (n < 0)
		throw invalid_argument(NegativeBitPosition);

	if (sign < 0) return (!Not().TestBit(n));

	wordNum = n / 32;
	if (wordNum >= magnitude.size()) return false;

	word = magnitude[magnitude.size() - 1 - wordNum];

	return ((uint32_t(word) >> (n & 31)) & 1) > 0;
}

BigInteger BigInteger::SetBit(const int32_t n) const
{
	if (n < 0)
		throw invalid_argument(InvalidBitAddress);

	if (TestBit(n)) return *this;

	BigInteger temp = *this;
	// TODO Handle negative values and zero
	if ((sign > 0) && (n < (temp.GetBitLength() - 1)))
		return FlipExistingBit(n);

	return Or(One.ShiftLeft(n));
}

BigInteger BigInteger::ClearBit(const int32_t n) const
{
	if (n < 0)
		throw invalid_argument(InvalidBitAddress);

	if (!TestBit(n)) return *this;

	BigInteger temp = *this;
	// TODO Handle negative values
	if ((sign > 0) && (n < (temp.GetBitLength() - 1)))
		return FlipExistingBit(n);

	return AndNot(One.ShiftLeft(n));
}

BigInteger BigInteger::FlipBit(const int32_t n) const
{
	if (n < 0)
		throw invalid_argument(InvalidBitAddress);

	BigInteger temp = *this;
	// TODO Handle negative values and zero
	if ((sign > 0) && (n < (temp.GetBitLength() - 1)))
		return FlipExistingBit(n);

	return Xor(One.ShiftLeft(n));
}

int32_t BigInteger::GetLowestSetBit() const
{
	if (sign == 0) return -1;

	return GetLowestSetBitMaskFirst(-1);
}

int ui64toa_s(uint64_t value, char* buffer, size_t size, int radix)
{
	if (!buffer || size == 0 || radix < 2 || radix > 36) {
		return EINVAL; // Invalid argument 
	}
	std::string result;
	do
	{
		int digit = value % radix;
		value /= radix;
		result.push_back(digit < 10 ? '0' + digit : 'a' + digit - 10);
	} while (value > 0);
	if (result.size() + 1 > size)
	{
		return ERANGE; // Buffer too small
	}
	std::reverse(result.begin(), result.end());
	std::copy(result.begin(), result.end(), buffer);
	buffer[result.size()] = '\0';
	return 0; // Success
}

void BigInteger::ToString(ostringstream &sl, int32_t radix, vector<BigInteger> &moduli, int32_t scale, BigInteger &pos)
{
	char temp[25];
	vector<BigInteger> qr;

	if (pos.GetBitLength() < 64)
	{
		ui64toa_s(uint64_t(pos.GetInt64Value()), &(temp[0]), 25, radix);

		uint32_t len = sl.str().size();
		if ((len > 1) || ((len == 1) && (sl.str()[0] != '-')))
		{
			AppendZeroExtendedString(sl, temp, 1 << scale);
		}
		else if (pos.GetSignValue() != 0)
		{
			sl << temp;
		}

		return;
	}

	--scale;
	qr = pos.DivideAndRemainder(moduli[scale]);

	ToString(sl, radix, moduli, scale, qr[0]);
	ToString(sl, radix, moduli, scale, qr[1]);
}

string BigInteger::ToString() const
{
	return ToString(10);
}

string BigInteger::ToString(const int32_t radix) const
{
	int32_t firstNonZero, pos, mask, bits, i, scale;
	ostringstream sl;
	vector<string> s;
	BigInteger u, q, r;
	vector<BigInteger> moduli;

	// TODO Make this method work for other radices (ideally 2 <= radix <= 36 as in Java)
	switch (radix)
	{
	case 2:
	case 8:
	case 10:
	case 16:
		// do nothing because it is in valid supported range
		break;

	default:
		throw invalid_argument(UnSupportedBase);
	}

	// NB: Can only happen to internally managed instances
	if ((!IsInitialized) && (magnitude.empty())) return "Null";

	if (sign == 0) return "0";

	// NOTE: This *should* be unnecessary, since the magnitude *should* never have leading zero digits
	firstNonZero = 0;
	while (firstNonZero < magnitude.size())
	{
		if (magnitude[firstNonZero] != 0) break;

		firstNonZero++;
	}

	if (firstNonZero == magnitude.size()) return "0";

	if (sign == -1)
		sl << '-';

	switch (radix)
	{
	case 2:
		pos = firstNonZero;

		sl << BigInteger::IntToBin(int32_t(magnitude[pos]));
		pos++;
		while (pos < magnitude.size())
		{
			AppendZeroExtendedString(sl, BigInteger::IntToBin(int64_t(magnitude[pos])), 32);
			pos++;
		}
		break;

	case 8:
		mask = (1 << 30) - 1;
		u = Abs();
		bits = u.GetBitLength();

		while (bits > 30)
		{
			s.push_back(BigInteger::IntToOctal(int32_t(u.GetInt32Value() & mask)));
			u = u.ShiftRight(30);
			bits = bits - 30;
		}

		sl << BigInteger::IntToOctal(int32_t(u.GetInt32Value()));
		i = s.size() - 1;
		while (i >= 0)
		{
			AppendZeroExtendedString(sl, s[i], 10);
			i--;
		}
		break;

	case 16:
		char temp[25];
		pos = firstNonZero;

		ui64toa_s(uint64_t(magnitude[pos]), &(temp[0]), 25, radix);
		sl << temp;
		pos++;

		while (pos < magnitude.size())
		{
			ui64toa_s(uint64_t(magnitude[pos]), &(temp[0]), 25, radix);
			AppendZeroExtendedString(sl, &(temp[0]), 8);
			pos++;
		}
		break;

		// TODO This could work for other radices if there is an alternative to Convert.ToString method
		// default:
	case 10:
		q = Abs();
		if (q.GetBitLength() < 64)
		{
			sl << q.GetInt64Value();
			return sl.str();
		}

		// TODO possibly cache power/exponent against radix?
		r = ValueOf(radix);
		while (r.CompareTo(q) <= 0)
		{
			moduli.push_back(r);
			r = r.Square();
		}

		scale = moduli.size();

		ToString(sl, radix, moduli, scale, q);
	}

	return sl.str();
}

bool BigInteger::Equals(const BigInteger &other) const
{
	return (sign == other.sign) && IsEqualMagnitude(other);
}

int32_t BigInteger::GetHashCode() const
{
	int32_t hc;

	hc = magnitude.size();
	if (magnitude.size() > 0)
	{
		hc = hc ^ magnitude[0];

		if (magnitude.size() > 1)
			hc = hc ^ magnitude[magnitude.size() - 1];
	}

	if (sign < 0) return ~hc;

	return hc;
}

int32_t BigInteger::BitCnt(const int32_t i)
{
	uint32_t u = uint32_t(i);

	u = u - ((u >> 1) & 0x55555555);
	u = (u & 0x33333333) + ((u >> 2) & 0x33333333);
	u = (u + (u >> 4)) & 0x0F0F0F0F;
	u = u + (u >> 8);
	u = u + (u >> 16);
	u = u & 0x3F;

	return int32_t(u);
}

int32_t BigInteger::BitLen(const int32_t w)
{
	uint32_t v, t;

	v = uint32_t(w);
	t = v >> 24;
	if (t != 0) return 24 + BitLengthTable[t];

	t = v >> 16;
	if (t != 0) return 16 + BitLengthTable[t];

	t = v >> 8;
	if (t != 0) return 8 + BitLengthTable[t];

	return BitLengthTable[v];
}

BigInteger BigInteger::ProbablePrime(const int32_t BitLength, IRandom random)
{
	return BigInteger(BitLength, 100, random);
}

BigInteger BigInteger::ValueOf(const int64_t value)
{
	if ((value >= 0) && (value < SMALL_CONSTANTS.size()))
		return SMALL_CONSTANTS[value];

	return CreateValueOf(value);
}

BigInteger BigInteger::Arbitrary(const int32_t sizeInBits)
{
	return BigInteger(sizeInBits, RandomSource);
}

void BigInteger::Boot()
{
	uint32_t i;
	vector<int32_t> primeList;
	int32_t product, j;
	BigInteger b;

	ZeroEncoding = vector<uint8_t>(0);
	ZeroMagnitude = vector<int32_t>(0);

	primeLists = vector<vector<int32_t>>({
		vector<int32_t>({ 3, 5, 7, 11, 13, 17, 19, 23 }),
		vector<int32_t>({ 29, 31, 37, 41, 43 }),
		vector<int32_t>({ 47, 53, 59, 61, 67 }),
		vector<int32_t>({ 71, 73, 79, 83 }),
		vector<int32_t>({ 89, 97, 101, 103 }),

		vector<int32_t>({ 107, 109, 113, 127 }),
		vector<int32_t>({ 131, 137, 139, 149 }),
		vector<int32_t>({ 151, 157, 163, 167 }),
		vector<int32_t>({ 173, 179, 181, 191 }),
		vector<int32_t>({ 193, 197, 199, 211 }),

		vector<int32_t>({ 223, 227, 229 }),
		vector<int32_t>({ 233, 239, 241 }),
		vector<int32_t>({ 251, 257, 263 }),
		vector<int32_t>({ 269, 271, 277 }),
		vector<int32_t>({ 281, 283, 293 }),

		vector<int32_t>({ 307, 311, 313 }),
		vector<int32_t>({ 317, 331, 337 }),
		vector<int32_t>({ 347, 349, 353 }),
		vector<int32_t>({ 359, 367, 373 }),
		vector<int32_t>({ 379, 383, 389 }),

		vector<int32_t>({ 397, 401, 409 }),
		vector<int32_t>({ 419, 421, 431 }),
		vector<int32_t>({ 433, 439, 443 }),
		vector<int32_t>({ 449, 457, 461 }),
		vector<int32_t>({ 463, 467, 479 }),

		vector<int32_t>({ 487, 491, 499 }),
		vector<int32_t>({ 503, 509, 521 }),
		vector<int32_t>({ 523, 541, 547 }),
		vector<int32_t>({ 557, 563, 569 }),
		vector<int32_t>({ 571, 577, 587 }),

		vector<int32_t>({ 593, 599, 601 }),
		vector<int32_t>({ 607, 613, 617 }),
		vector<int32_t>({ 619, 631, 641 }),
		vector<int32_t>({ 643, 647, 653 }),
		vector<int32_t>({ 659, 661, 673 }),

		vector<int32_t>({ 677, 683, 691 }),
		vector<int32_t>({ 701, 709, 719 }),
		vector<int32_t>({ 727, 733, 739 }),
		vector<int32_t>({ 743, 751, 757 }),
		vector<int32_t>({ 761, 769, 773 }),

		vector<int32_t>({ 787, 797, 809 }),
		vector<int32_t>({ 811, 821, 823 }),
		vector<int32_t>({ 827, 829, 839 }),
		vector<int32_t>({ 853, 857, 859 }),
		vector<int32_t>({ 863, 877, 881 }),

		vector<int32_t>({ 883, 887, 907 }),
		vector<int32_t>({ 911, 919, 929 }),
		vector<int32_t>({ 937, 941, 947 }),
		vector<int32_t>({ 953, 967, 971 }),
		vector<int32_t>({ 977, 983, 991 }),

		vector<int32_t>({ 997, 1009, 1013 }),
		vector<int32_t>({ 1019, 1021, 1031 }),
		vector<int32_t>({ 1033, 1039, 1049 }),
		vector<int32_t>({ 1051, 1061, 1063 }),
		vector<int32_t>({ 1069, 1087, 1091 }),

		vector<int32_t>({ 1093, 1097, 1103 }),
		vector<int32_t>({ 1109, 1117, 1123 }),
		vector<int32_t>({ 1129, 1151, 1153 }),
		vector<int32_t>({ 1163, 1171, 1181 }),
		vector<int32_t>({ 1187, 1193, 1201 }),

		vector<int32_t>({ 1213, 1217, 1223 }),
		vector<int32_t>({ 1229, 1231, 1237 }),
		vector<int32_t>({ 1249, 1259, 1277 }),
		vector<int32_t>({ 1279, 1283, 1289 })
	});

	SecureRandom::Boot();

	RandomSource = make_shared<SecureRandom>();

	Zero = BigInteger(0, ZeroMagnitude, false);
	Zero.nBits = 0;
	Zero.nBitLength = 0;

	SMALL_CONSTANTS = vector<BigInteger>(17);

	SMALL_CONSTANTS[0] = Zero;

	i = 1;
	while (i < uint32_t(SMALL_CONSTANTS.size()))
	{
		SMALL_CONSTANTS[i] = CreateUValueOf(i);
		i++;
	}

	One = SMALL_CONSTANTS[1];
	Two = SMALL_CONSTANTS[2];
	Three = SMALL_CONSTANTS[3];
	Ten = SMALL_CONSTANTS[10];

	radix2 = ValueOf(2);
	radix2E = radix2.Pow(chunk2);

	radix8 = ValueOf(8);
	radix8E = radix8.Pow(chunk8);

	radix10 = ValueOf(10);

	radix10E = radix10.Pow(chunk10);

	radix16 = ValueOf(16);
	radix16E = radix16.Pow(chunk16);

	primeProducts = vector<int32_t>(primeLists.size());

	for (i = 0; i < primeLists.size(); i++)
	{
		primeList = primeLists[i];
		product = primeList[0];
		for (j = 1; j < primeList.size(); j++)
			product = product * primeList[j];

		primeProducts[i] = product;
	}

}

int32_t BigInteger::GetBitLength()
{
	if (sign == 0) nBitLength = 0;
	else
		nBitLength = CalcBitLength(sign, 0, magnitude);

	return nBitLength;
}

int32_t BigInteger::GetBitCount()
{
	int32_t sum, i;

	if (sign < 0)
		// TODO Optimise this case
		nBits = Not().GetBitCount();
	else
	{
		sum = 0;
		for (i = 0; i < magnitude.size(); i++)
			sum = sum + BitCnt(magnitude[i]);

		nBits = sum;
	}

	return nBits;
}

int32_t BigInteger::GetInt32Value() const
{
	int32_t n, v;

	if (sign == 0) return 0;

	n = magnitude.size();

	v = magnitude[n - 1];

	if (sign < 0) return -v;

	return v;
}

int64_t BigInteger::GetInt64Value() const
{
	int32_t n;
	int64_t v;

	if (sign == 0) return 0;

	n = magnitude.size();

	v = magnitude[n - 1] & IMASK;
	if (n > 1)
		v = v | ((magnitude[n - 2] & IMASK) << 32);

	if (sign < 0) return -v;

	return v;
}


BigInteger BigInteger::AddToMagnitude(const vector<int32_t> &magToAdd) const
{
	vector<int32_t> big, small, bigCopy;
	uint32_t limit;
	bool possibleOverflow;

	if (magnitude.size() < magToAdd.size())
	{
		big = magToAdd;
		small = magnitude;
	}
	else
	{
		big = magnitude;
		small = magToAdd;
	}

	// Conservatively avoid over-allocation when no overflow possible
	limit = UINT32_MAX;
	if (big.size() == small.size())
		limit = limit - uint32_t(small[0]);

	possibleOverflow = uint32_t(big[0]) >= limit;

	if (possibleOverflow)
	{
		bigCopy = vector<int32_t>(big.size() + 1);
		memmove(&bigCopy[1], &big[0], big.size() * sizeof(int32_t));
	}
	else
		bigCopy = big;

	bigCopy = AddMagnitudes(bigCopy, small);

	return BigInteger(sign, bigCopy, possibleOverflow);
}

bool BigInteger::QuickPow2Check() const
{
	BigInteger temp = *this;
	return (sign > 0) && (temp.GetBitCount() == 1);
}

vector<int32_t> BigInteger::Divide(vector<int32_t> &x, vector<int32_t> &y) const
{
	int32_t xStart, yStart, xyCmp, yBitLength, xBitLength, shift,
		iCountStart, cBitLength, cStart, len;
	vector<int32_t> Count, iCount, c;
	uint32_t firstC, firstX;

	xStart = 0;
	while ((xStart < x.size()) && (x[xStart] == 0)) xStart++;

	yStart = 0;
	while ((yStart < y.size()) && (y[yStart] == 0)) yStart++;

	xyCmp = CompareNoLeadingZeroes(xStart, x, yStart, y);

	if (xyCmp > 0)
	{
		yBitLength = CalcBitLength(1, yStart, y);
		xBitLength = CalcBitLength(1, xStart, x);
		shift = xBitLength - yBitLength;

		iCountStart = 0;

		cStart = 0;
		cBitLength = yBitLength;
		if (shift > 0)
		{
			iCount = vector<int32_t>((shift >> 5) + 1);

			iCount[0] = 1 << (shift & 31);

			c = ShiftLeft(y, shift);
			cBitLength = cBitLength + shift;
		}
		else
		{
			iCount = vector<int32_t>({ 1 });

			len = y.size() - yStart;
			c = vector<int32_t>(len);
			memmove(&c[0], &y[yStart], len * sizeof(int32_t));
		}

		Count = vector<int32_t>(iCount.size());

		while (true)
		{
			if ((cBitLength < xBitLength) || (CompareNoLeadingZeroes(xStart, x, cStart, c) >= 0))
			{
				Subtract(xStart, x, cStart, c);
				AddMagnitudes(Count, iCount);

				while (x[xStart] == 0)
				{
					xStart++;
					if (xStart == x.size()) return Count;
				}

				xBitLength = 32 * (x.size() - xStart - 1) + BitLen(x[xStart]);

				if (xBitLength <= yBitLength)
				{
					if (xBitLength < yBitLength) return Count;

					xyCmp = CompareNoLeadingZeroes(xStart, x, yStart, y);

					if (xyCmp <= 0) break;
				}
			}

			shift = cBitLength - xBitLength;

			// NB: The case where c[cStart] is 1-bit is harmless
			if (shift == 1)
			{
				firstC = uint32_t(uint32_t(c[cStart]) >> 1);
				firstX = uint32_t(x[xStart]);
				if (firstC > firstX) shift++;
			}

			if (shift < 2)
			{
				ShiftRightOneInPlace(cStart, c);
				cBitLength--;
				ShiftRightOneInPlace(iCountStart, iCount);
			}
			else
			{
				ShiftRightInPlace(cStart, c, shift);
				cBitLength = cBitLength - shift;
				ShiftRightInPlace(iCountStart, iCount, shift);
			}

			while (c[cStart] == 0) cStart++;

			while (iCount[iCountStart] == 0) iCountStart++;
		}
	}
	else
		Count = vector<int32_t>(1);

	if (xyCmp == 0)
	{
		AddMagnitudes(Count, One.magnitude);
		memset(&x[xStart], 0, (x.size() - xStart) * sizeof(int32_t));
	}

	return Count;
}

bool BigInteger::IsEqualMagnitude(const BigInteger &x) const
{
	vector<int32_t> xMag;
	int i;

	xMag = x.magnitude;
	if (magnitude.size() != x.magnitude.size()) return false;

	for (i = 0; i < magnitude.size(); i++)
	{
		if (magnitude[i] != x.magnitude[i]) return false;
	}

	return true;
}

bool BigInteger::CheckProbablePrime(const int32_t certainty, IRandom random, const bool randomlySelected)
{
	int32_t numLists, i, j, prime, qRem, test;
	vector<int32_t> primeList;

	// Try to reduce the penalty for really small numbers
	numLists = std::min(GetBitLength() - 1, (int)primeLists.size());

	for (i = 0; i < numLists; i++)
	{
		test = Remainder(primeProducts[i]);

		primeList = primeLists[i];

		for (j = 0; j < primeList.size(); j++)
		{
			prime = primeList[j];
			qRem = test % prime;
			if (qRem == 0)
				// We may find small numbers in the list
				return (GetBitLength() < 16) && (GetInt32Value() == prime);
		}

	}

	// TODO Is it worth trying to create a hybrid of these two?
	return RabinMillerTest(certainty, random, randomlySelected);
}

BigInteger BigInteger::ModInversePow2(const BigInteger &value) const
{
	int32_t Pow, bitsCorrect;
	int64_t inv64;
	BigInteger x, d, t, m = value;

	if (!TestBit(0))
		throw invalid_argument(NotRelativelyPrime);

	Pow = m.GetBitLength() - 1;

	inv64 = ModInverse64(GetInt64Value());
	if (Pow < 64) inv64 = inv64 & ((int64_t(1) << Pow) - 1);

	x = BigInteger::ValueOf(inv64);

	if (Pow > 64)
	{
		d = Remainder(m);
		bitsCorrect = 64;

		do
		{
			t = x.Multiply(d).Remainder(m);
			x = x.Multiply(Two.Subtract(t)).Remainder(m);
			bitsCorrect = bitsCorrect << 1;
		} while (bitsCorrect < Pow);
	}

	if (x.sign < 0) x = x.Add(m);

	return x;
}

int32_t BigInteger::GetMQuote()
{
	int32_t d;

	d = int32_t(magnitude[magnitude.size() - 1]);
	d = -d;

	mQuote = ModInverse32(d);

	return mQuote;
}

vector<int32_t> BigInteger::Remainder(vector<int32_t> &x, const vector<int32_t> &y) const
{
	int32_t xStart, yStart, xyCmp, yBitLength, xBitLength,
		shift, cBitLength, cStart, len;
	vector<int32_t> c;
	uint32_t firstC, firstX;

	xStart = 0;
	while ((xStart < x.size()) && (x[xStart] == 0)) xStart++;

	yStart = 0;
	while ((yStart < y.size()) && (y[yStart] == 0)) yStart++;

	xyCmp = CompareNoLeadingZeroes(xStart, x, yStart, y);

	if (xyCmp > 0)
	{
		yBitLength = CalcBitLength(1, yStart, y);
		xBitLength = CalcBitLength(1, xStart, x);
		shift = xBitLength - yBitLength;

		cStart = 0;
		cBitLength = yBitLength;
		if (shift > 0)
		{
			c = ShiftLeft(y, shift);
			cBitLength = cBitLength + shift;
		}
		else
		{
			len = y.size() - yStart;
			c = vector<int32_t>(len);
			memmove(&c[0], &y[yStart], len * sizeof(int32_t));
		}

		while (true)
		{
			if ((cBitLength < xBitLength) || (CompareNoLeadingZeroes(xStart, x, cStart, c) >= 0))
			{
				Subtract(xStart, x, cStart, c);

				while (x[xStart] == 0)
				{
					xStart++;
					if (xStart == x.size()) return x;
				}

				xBitLength = 32 * (x.size() - xStart - 1) + BitLen(x[xStart]);

				if (xBitLength <= yBitLength)
				{
					if (xBitLength < yBitLength) return x;

					xyCmp = CompareNoLeadingZeroes(xStart, x, yStart, y);

					if (xyCmp <= 0) break;
				}

			}

			shift = cBitLength - xBitLength;

			// NB: The case where c[cStart] is 1-bit is harmless
			if (shift == 1)
			{
				firstC = uint32_t(uint32_t(c[cStart]) >> 1);
				firstX = uint32_t(x[xStart]);
				if (firstC > firstX) shift++;
			}

			if (shift < 2)
			{
				ShiftRightOneInPlace(cStart, c);
				cBitLength--;
			}
			else
			{
				ShiftRightInPlace(cStart, c, shift);
				cBitLength = cBitLength - shift;
			}

			while (c[cStart] == 0) cStart++;
		}
	}

	if (xyCmp == 0)
		memset(&x[xStart], 0, (x.size() - xStart) * sizeof(int32_t));

	return x;
}

vector<int32_t> BigInteger::LastNBits(const int32_t n) const
{
	int32_t numWords, excessBits;
	vector<int32_t> result;

	if (n < 1) return ZeroMagnitude;

	numWords = (n + BitsPerInt - 1) / BitsPerInt;
	numWords = std::min(numWords, (int)magnitude.size());

	result = vector<int32_t>(numWords);
	memmove(&result[0], &magnitude[magnitude.size() - numWords], numWords * sizeof(int32_t));

	excessBits = (numWords << 5) - n;
	if (excessBits > 0)
		result[0] = result[0] & (int32_t(uint32_t(UINT32_MAX) >> excessBits));

	return result;
}

BigInteger BigInteger::DivideWords(const int32_t w) const
{
	int32_t n;
	vector<int32_t> mag;

	n = magnitude.size();
	if (w >= n) return Zero;

	mag = vector<int32_t>(n - w);
	memmove(&mag[0], &magnitude[0], (n - w) * sizeof(int32_t));

	return BigInteger(sign, mag, false);
}

BigInteger BigInteger::RemainderWords(const int32_t w) const
{
	int32_t n;
	vector<int32_t> mag;

	n = magnitude.size();
	if (w >= n) return *this;

	mag = vector<int32_t>(w);
	memmove(&mag[0], &magnitude[n - w], w * sizeof(int32_t));

	return BigInteger(sign, mag, false);
}

vector<uint8_t> BigInteger::ToByteArray(const bool _unsigned) const
{
	int32_t nBits, nBytes, magIndex, bytesIndex;
	uint32_t mag, lastMag;
	vector<uint8_t> bytes;
	bool carry;

	if (sign == 0)
	{
		if (_unsigned) return ZeroEncoding;

		return vector<uint8_t>(1);
	}

	BigInteger temp = *this;
	if ((_unsigned) && (sign > 0)) nBits = temp.GetBitLength();
	else
		nBits = temp.GetBitLength() + 1;

	nBytes = GetByteLength(nBits);
	bytes = vector<uint8_t>(nBytes);

	magIndex = magnitude.size();
	bytesIndex = bytes.size();

	if (sign > 0)
	{
		while (magIndex > 1)
		{
			magIndex--;
			mag = uint32_t(magnitude[magIndex]);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 8);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 16);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 24);
		}

		lastMag = uint32_t(magnitude[0]);
		while (lastMag > UINT8_MAX)
		{
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(lastMag);
			lastMag = lastMag >> 8;
		}

		bytesIndex--;
		bytes[bytesIndex] = uint8_t(lastMag);
	}
	else // sign < 0
	{
		carry = true;

		while (magIndex > 1)
		{
			magIndex--;
			mag = ~(uint32_t(magnitude[magIndex]));

			if (carry)
			{
				mag++;
				carry = (mag == INT32_MIN);
			}

			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 8);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 16);
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(mag >> 24);
		}

		lastMag = uint32_t(magnitude[0]);

		if (carry)
			// Never wraps because magnitude[0] != 0
			lastMag--;

		while (lastMag > UINT8_MAX)
		{
			bytesIndex--;
			bytes[bytesIndex] = uint8_t(~lastMag);
			lastMag = lastMag >> 8;
		}

		bytesIndex--;
		bytes[bytesIndex] = uint8_t(~lastMag);

		if (bytesIndex > 0)
		{
			bytesIndex--;
			bytes[bytesIndex] = UINT8_MAX;
		}

	}

	return bytes;
}

BigInteger BigInteger::FlipExistingBit(const int32_t n) const
{
	vector<int32_t> mag;

	mag = magnitude;
	mag[mag.size() - 1 - (uint32_t(n) >> 5)] = mag[mag.size() - 1 - (uint32_t(n) >> 5)] ^ (1 << (n & 31));

	return BigInteger(sign, mag, false);
}

int32_t BigInteger::GetLowestSetBitMaskFirst(const int32_t firstWordMask) const
{
	int32_t w, offset;
	uint32_t word;

	w = magnitude.size();
	offset = 0;

	w--;
	word = uint32_t(magnitude[w] & firstWordMask);

	while (word == 0)
	{
		w--;
		word = uint32_t(magnitude[w]);
		offset = offset + 32;
	}

	while ((word & 0xFF) == 0)
	{
		word = word >> 8;
		offset = offset + 8;
	}

	while ((word & 1) == 0)
	{
		word = word >> 1;
		offset++;
	}

	return offset;
}

void BigInteger::ParseString(const string &str, const int32_t radix)
{
	auto style = 0;
	int32_t chunk, index, Next;
	BigInteger r, rE, b, bi;
	string dVal, s, temp;
	uint64_t i;

	if (str.size() == 0)
		throw invalid_argument(ZeroLengthBigInteger);

	sign = -1;
	nBits = -1;
	nBitLength = -1;
	mQuote = 0;
	IsInitialized = true;

	switch (radix)
	{
	case 2:
		// Is there anyway to restrict to binary digits?
		style = NumberStyles::Integer;
		chunk = chunk2;
		r = radix2;
		rE = radix2E;
		break;

	case 8:
		// Is there anyway to restrict to octal digits?
		style = NumberStyles::Integer;
		chunk = chunk8;
		r = radix8;
		rE = radix8E;
		break;

	case 10:
		// This style seems to handle spaces and minus sign already (our processing redundant?)
		style = NumberStyles::Integer;
		chunk = chunk10;
		r = radix10;
		rE = radix10E;
		break;

	case 16:
		// TODO Should this be HexNumber?
		style = NumberStyles::AllowHexSpecifier;
		chunk = chunk16;
		r = radix16;
		rE = radix16E;
		break;

	default:
		throw invalid_argument(InvalidBase);
	}

	index = 0;
	sign = 1;

	if (str[0] == '-')
	{
		if (str.size() == 1)
			throw invalid_argument(ZeroLengthBigInteger);

		sign = -1;
		index = 1;
	}

	// strip leading zeros from the string str
	while (index < str.size())
	{
		dVal = str[index];

		temp = dVal;

		if (strtoull(temp.c_str(), NULL, radix) == 0) index++;
		else
			break;
	}

	if (index >= str.size())
	{
		// zero value - we're done
		sign = 0;
		magnitude = ZeroMagnitude;
		return;
	}

	/// ///
	// could we work out the max number of ints required to store
	// str.Length digits in the given base, then allocate that
	// storage in one hit?, then Generate the magnitude in one hit too?
	/// ///

	b = Zero;

	Next = index + chunk;

	while (Next <= str.size())
	{
		s = str.substr(index, chunk);

		temp = s;

		i = strtoull(temp.c_str(), NULL, radix);

		bi = CreateUValueOf(i);

		switch (radix)
		{
		case 2:
			// TODO Need this because we are parsing in radix 10 above
			if (i >= 2)
				throw invalid_argument(BadCharacterRadix2);

			// TODO Parse 64 bits at a time
			b = b.ShiftLeft(1);
			break;

		case 8:
			// TODO Need this because we are parsing in radix 10 above
			if (i >= 8)
				throw invalid_argument(BadCharacterRadix8);

			// TODO Parse 63 bits at a time
			b = b.ShiftLeft(3);
			break;

		case 16:
			b = b.ShiftLeft(64);
			break;

		default:
			b = b.Multiply(rE);
		}

		b = b.Add(bi);

		index = Next;
		Next = Next + chunk;
	}

	if (index < str.size())
	{
		s = str.substr(index, str.size() - (index - 1));

		temp = s;

		i = strtoull(temp.c_str(), NULL, radix);

		bi = CreateUValueOf(i);

		if (b.sign > 0)
		{
			if (radix == 2)
				// TODO Parse all bits at once
				b = b.ShiftLeft(s.size());
			else if (radix == 8)
				// TODO Parse all bits at once
				b = b.ShiftLeft(s.size() * 3);
			else if (radix == 16) b = b.ShiftLeft(s.size() << 2);
			else
				b = b.Multiply(r.Pow(s.size()));

			b = b.Add(bi);
		}
		else
			b = bi;
	}

	magnitude = b.magnitude;
}

void BigInteger::ParseBytes(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length)
{
	int32_t endPoint, iBval, numBytes, index;
	vector<uint8_t> inverse;

	if (length == 0)
		throw invalid_argument(ZeroLengthBigInteger);

	sign = -1;
	nBits = -1;
	nBitLength = -1;
	mQuote = 0;
	IsInitialized = true;

	// TODO Move this processing into MakeMagnitude (provide sign argument)
	if (int8_t(bytes[offset]) < 0)
	{
		sign = -1;

		endPoint = offset + length;

		iBval = offset;

		// strip leading sign bytes
		while ((iBval < endPoint) && (int8_t(bytes[iBval]) == -1)) iBval++;

		if (iBval >= endPoint) magnitude = One.magnitude;
		else
		{
			numBytes = endPoint - iBval;
			inverse = vector<uint8_t>(numBytes);

			index = 0;
			while (index < numBytes)
			{
				inverse[index] = uint8_t(~bytes[iBval]);
				index++;
				iBval++;
			}

			index--;
			while (inverse[index] == UINT8_MAX)
			{
				inverse[index] = INT8_MIN;
				index--;
			}

			inverse[index] = uint8_t(inverse[index] + 1);

			magnitude = MakeMagnitude(inverse, 0, inverse.size());
		}
	}
	else
	{
		// strip leading zero bytes and return magnitude bytes
		magnitude = MakeMagnitude(bytes, offset, length);

		if (magnitude.size() > 0) sign = 1;
		else
			sign = 0;
	}

}

void BigInteger::ParseBytesWithSign(const int32_t _sign, const vector<uint8_t> &bytes, const int32_t offset, const int32_t length)
{
	if ((_sign < -1) || (_sign > 1))
		throw invalid_argument(InvalidSign);

	sign = -1;
	nBits = -1;
	nBitLength = -1;
	mQuote = 0;
	IsInitialized = true;

	if (_sign == 0)
	{
		sign = 0;
		magnitude = ZeroMagnitude;
	}
	else
	{
		// copy bytes
		magnitude = MakeMagnitude(bytes, offset, length);
		if (magnitude.size() < 1) sign = 0;
		else
			sign = _sign;
	}

}

int32_t BigInteger::GetByteLength(const int32_t nBits)
{
	return (nBits + BitsPerByte - 1) / BitsPerByte;
}

vector<int32_t> BigInteger::MakeMagnitude(const vector<uint8_t> &bytes, const int32_t offset, const int32_t length)
{
	int32_t endPoint, firstSignificant, nInts, bCount, v, magnitudeIndex, i;
	vector<int32_t> mag;

	endPoint = offset + length;

	// strip leading zeros
	firstSignificant = offset;
	while ((firstSignificant < endPoint) && (bytes[firstSignificant] == 0))
		firstSignificant++;

	if (firstSignificant >= endPoint) return ZeroMagnitude;

	nInts = (endPoint - firstSignificant + 3) / BytesPerInt;
	bCount = (endPoint - firstSignificant) % BytesPerInt;
	if (bCount == 0) bCount = BytesPerInt;

	if (nInts < 1) return ZeroMagnitude;

	mag = vector<int32_t>(nInts);

	v = 0;
	magnitudeIndex = 0;

	i = firstSignificant;
	while (i < endPoint)
	{
		v = v << 8;
		v = v | (bytes[i] & 0xFF);
		bCount--;
		if (bCount <= 0)
		{
			mag[magnitudeIndex] = v;
			magnitudeIndex++;
			bCount = BytesPerInt;
			v = 0;
		}

		i++;
	}

	if (magnitudeIndex < mag.size()) mag[magnitudeIndex] = v;

	return mag;
}

vector<int32_t> BigInteger::AddMagnitudes(vector<int32_t> &a, const vector<int32_t> &b)
{
	int32_t tI, vI;
	int64_t m;

	tI = a.size() - 1;
	vI = b.size() - 1;
	m = 0;

	while (vI >= 0)
	{
		m = m + (int64_t(uint32_t(a[tI])) + int64_t(uint32_t(b[vI])));
		vI--;
		a[tI] = int32_t(m);
		tI--;
		m = int64_t(uint64_t(uint64_t(m) >> 32));
	}

	if (m != 0)
	{
		while (tI >= 0)
		{
			a[tI] = a[tI] + 1;

			if (a[tI] != 0)	break;

			tI--;
		}

	}

	return a;
}

int32_t BigInteger::CalcBitLength(const int32_t sign, const int32_t _indx, const vector<int32_t> &mag)
{
	int32_t BitLength, firstMag, indx = _indx;

	while (true)
	{
		if (indx >= mag.size()) return 0;

		if (mag[indx] != 0) break;

		indx++;
	}

	// bit length for everything after the first int
	BitLength = 32 * ((mag.size() - indx) - 1);

	// and determine bitlength of first int
	firstMag = mag[indx];
	BitLength = BitLength + BitLen(firstMag);

	// Check for negative powers of two
	if ((sign < 0) && ((firstMag & int32_t(-firstMag)) == firstMag))
	{
		do
		{
			indx++;
			if (indx >= mag.size())
			{
				BitLength--;
				break;
			}
		} while (mag[indx] == 0);
	}

	return BitLength;
}

int32_t BigInteger::CompareTo(const int32_t _xIndx, const vector<int32_t> &x, const int32_t _yIndx, const vector<int32_t> &y)
{
	int32_t xIndx = _xIndx, yIndx = _yIndx;

	while ((xIndx != x.size()) && (x[xIndx] == 0)) xIndx++;

	while ((yIndx != y.size()) && (y[yIndx] == 0)) yIndx++;

	return CompareNoLeadingZeroes(xIndx, x, yIndx, y);
}

int32_t BigInteger::CompareNoLeadingZeroes(const int32_t _xIndx, const vector<int32_t> &x, const int32_t _yIndx, const vector<int32_t> &y)
{
	int32_t diff;
	uint32_t v1, v2, xIndx = _xIndx, yIndx = _yIndx;

	diff = (x.size() - y.size()) - (xIndx - yIndx);

	if (diff != 0)
	{
		if (diff < 0) return -1;

		return 1;
	}

	// lengths of magnitudes the same, test the magnitude values

	while (xIndx < x.size())
	{
		v1 = uint32_t(x[xIndx]);
		xIndx++;
		v2 = uint32_t(y[yIndx]);
		yIndx++;

		if (v1 != v2)
		{
			if (v1 < v2) return -1;

			return 1;
		}
	}

	return 0;
}

int32_t BigInteger::ModInverse32(const int32_t d)
{
	int32_t x;

	// Newton's method with initial estimate "correct to 4 bits"
	x = d + (((d + 1) & 4) << 1); // d.x == 1 mod 2**4
	x = x * (2 - (d * x)); // d.x == 1 mod 2**8
	x = x * (2 - (d * x)); // d.x == 1 mod 2**16
	x = x * (2 - (d * x)); // d.x == 1 mod 2**32

	return x;
}

int64_t BigInteger::ModInverse64(const int64_t d)
{
	int64_t x;

	// Newton's method with initial estimate "correct to 4 bits"
	x = d + (((d + int64_t(1)) & int64_t(4)) << 1); // d.x == 1 mod 2**4
	x = x * (2 - (d * x)); // d.x == 1 mod 2**8
	x = x * (2 - (d * x)); // d.x == 1 mod 2**16
	x = x * (2 - (d * x)); // d.x == 1 mod 2**32
	x = x * (2 - (d * x)); // d.x == 1 mod 2**64

	return x;
}

BigInteger BigInteger::ExtEuclid(const BigInteger &a, const BigInteger &b, BigInteger &u1Out)
{
	BigInteger u1, v1, u3, v3, oldU1;
	vector<BigInteger> q;

	u1 = One;
	v1 = Zero;
	u3 = a;
	v3 = b;

	if (v3.sign > 0)
	{
		while (true)
		{
			q = u3.DivideAndRemainder(v3);
			u3 = v3;
			v3 = q[1];

			oldU1 = u1;
			u1 = v1;

			if (v3.sign <= 0) break;

			v1 = oldU1.Subtract(v1.Multiply(q[0]));
		}
	}

	u1Out = u1;

	return u3;
}

BigInteger BigInteger::ModPowBarrett(const BigInteger &b, const BigInteger &_e, const BigInteger &_m)
{
	int32_t k, extraBits, expLength, numPowers, i, window, mult,
		lastZeroes, windowPos, bits, j;
	BigInteger mr, yu, b2, y, e = _e, m = _m;
	vector<BigInteger> oddPowers;
	vector<int32_t> windowList;

	k = m.magnitude.size();
	mr = One.ShiftLeft((k + 1) << 5);
	yu = One.ShiftLeft(k << 6).Divide(m);

	// Sliding window from MSW to LSW
	extraBits = 0;
	expLength = e.GetBitLength();
	while (expLength > ExpWindowThresholds[extraBits]) extraBits++;

	numPowers = 1 << extraBits;
	oddPowers = vector<BigInteger>(numPowers);
	oddPowers[0] = b;

	BigInteger bSquare = b.Square();
	b2 = ReduceBarrett(bSquare, m, mr, yu);

	for (i = 1; i < numPowers; i++)
	{

		BigInteger x = oddPowers[i - 1].Multiply(b2);
		oddPowers[i] = ReduceBarrett(x, m, mr, yu);
	}

	windowList = GetWindowList(e.magnitude, extraBits);
	window = windowList[0];
	mult = window & 0xFF;
	lastZeroes = window >> 8;

	if (mult == 1)
	{
		y = b2;
		lastZeroes--;
	}
	else
		y = oddPowers[mult >> 1];

	windowPos = 1;
	window = windowList[windowPos];
	windowPos++;
	while (window != -1)
	{
		mult = window & 0xFF;

		bits = lastZeroes + BitLengthTable[mult];

		j = 0;
		while (j < bits)
		{
			BigInteger ySquare = y.Square();
			y = ReduceBarrett(ySquare, m, mr, yu);
			j++;
		}


		BigInteger yMul = y.Multiply(oddPowers[mult >> 1]);
		y = ReduceBarrett(yMul, m, mr, yu);

		lastZeroes = window >> 8;

		window = windowList[windowPos];
		windowPos++;
	}

	i = 0;
	while (i < lastZeroes)
	{
		BigInteger ySq = y.Square();
		y = ReduceBarrett(ySq, m, mr, yu);
		i++;
	}

	return y;
}

BigInteger BigInteger::ReduceBarrett(BigInteger &x, BigInteger &m, const BigInteger &mr, const BigInteger &yu)
{
	int32_t xLen, mLen, k;
	BigInteger q1, q2, q3, r1, r2, r3;

	xLen = x.GetBitLength();
	mLen = m.GetBitLength();
	if (xLen < mLen) return x;

	if ((xLen - mLen) > 1)
	{
		k = m.magnitude.size();

		q1 = x.DivideWords(k - 1);
		q2 = q1.Multiply(yu); // TODO Only need partial multiplication here
		q3 = q2.DivideWords(k + 1);

		r1 = x.RemainderWords(k + 1);
		r2 = q3.Multiply(m); // TODO Only need partial multiplication here
		r3 = r2.RemainderWords(k + 1);

		x = r1.Subtract(r3);
		if (x.sign < 0) x = x.Add(mr);
	}

	while (x.CompareTo(m) >= 0)
		x = x.Subtract(m);

	return x;
}

BigInteger BigInteger::ModPowMonty(BigInteger &b, const BigInteger &_e, const BigInteger &_m, const bool convert)
{
	int32_t n, powR, extraBits, expLength, numPowers, i, window,
		mult, lastZeroes, windowPos, bits, j;
	bool smallMontyModulus;
	uint32_t mDash;
	vector<int32_t> yAccum, zVal, tmp, zSquared, windowList, yVal;
	vector<vector<int32_t>> oddPowers;
	BigInteger m = _m, e = _e;

	n = m.magnitude.size();
	powR = 32 * n;
	smallMontyModulus = (m.GetBitLength() + 2) <= powR;
	mDash = uint32_t(m.GetMQuote());

	// tmp = this * R mod m
	if (convert)
		b = b.ShiftLeft(powR).Remainder(m);

	yAccum = vector<int32_t>(n + 1);

	zVal = b.magnitude;
	if (zVal.size() < n)
	{
		tmp = vector<int32_t>(n);
		memmove(&tmp[n - zVal.size()], &zVal[0], zVal.size() * sizeof(int32_t));
		zVal = tmp;
	}

	// Sliding window from MSW to LSW

	extraBits = 0;

	// Filter the common case of small RSA exponents with few bits set
	if ((e.magnitude.size() > 1) || (e.GetBitCount() > 2))
	{
		expLength = e.GetBitLength();
		while (expLength > ExpWindowThresholds[extraBits])
			extraBits++;
	}

	numPowers = 1 << extraBits;

	oddPowers = vector<vector<int32_t>>(numPowers);
	oddPowers[0] = zVal;

	zSquared = zVal;

	SquareMonty(yAccum, zSquared, m.magnitude, mDash, smallMontyModulus);

	for (i = 1; i < numPowers; i++)
	{
		oddPowers[i] = oddPowers[i - 1];
		MultiplyMonty(yAccum, oddPowers[i], zSquared, m.magnitude, mDash, smallMontyModulus);

	}

	windowList = GetWindowList(e.magnitude, extraBits);
	window = windowList[0];
	mult = window & 0xFF;
	lastZeroes = window >> 8;

	if (mult == 1)
	{
		yVal = zSquared;
		lastZeroes--;
	}
	else
		yVal = oddPowers[mult >> 1];

	windowPos = 1;
	window = windowList[windowPos];
	windowPos++;
	while (int32_t(window) != int32_t(-1))
	{
		mult = window & 0xFF;

		bits = lastZeroes + BitLengthTable[mult];

		j = 0;
		while (j < bits)
		{
			SquareMonty(yAccum, yVal, m.magnitude, mDash, smallMontyModulus);
			j++;
		}

		MultiplyMonty(yAccum, yVal, oddPowers[mult >> 1], m.magnitude, mDash, smallMontyModulus);

		lastZeroes = window >> 8;
		window = windowList[windowPos];
		windowPos++;
	}

	i = 0;
	while (i < lastZeroes)
	{
		SquareMonty(yAccum, yVal, m.magnitude, mDash, smallMontyModulus);
		i++;
	}

	if (convert) MontgomeryReduce(yVal, m.magnitude, mDash);
	else if ((smallMontyModulus) && (CompareTo(0, yVal, 0, m.magnitude) >= 0))
		Subtract(0, yVal, 0, m.magnitude);

	return BigInteger(1, yVal, true);
}

vector<int32_t> BigInteger::GetWindowList(const vector<int32_t> &mag, const int32_t extraBits)
{
	int32_t i, v, leadingBits, resultSize, resultPos, bitPos, mult, multLimit, zeroes;

	v = mag[0];
	leadingBits = BitLen(v);

	resultSize = (((mag.size() - 1) << 5) + leadingBits) / (1 + extraBits) + 2;
	vector<int32_t> result = vector<int32_t>(resultSize);
	resultPos = 0;

	bitPos = 33 - leadingBits;
	v = v << bitPos;

	mult = 1;
	multLimit = 1 << extraBits;
	zeroes = 0;

	i = 0;
	while (true)
	{
		while (bitPos < 32)
		{
			if (mult < multLimit)
				mult = (mult << 1) | int32_t((uint32_t(v) >> 31));
			else if (v < 0)
			{
				result[resultPos] = CreateWindowEntry(mult, zeroes);
				resultPos++;
				mult = 1;
				zeroes = 0;
			}
			else
				zeroes++;

			v = v << 1;
			bitPos++;
		}

		i++;
		if (i == mag.size())
		{
			result[resultPos] = CreateWindowEntry(mult, zeroes);
			resultPos++;
			break;
		}

		v = mag[i];
		bitPos = 0;
	}

	result[resultPos] = -1;

	return result;
}

int32_t BigInteger::CreateWindowEntry(int32_t mult, int32_t zeroes)
{
	while ((mult & 1) == 0)
	{
		mult = mult >> 1;
		zeroes++;
	}

	return mult | (zeroes << 8);
}

vector<int32_t> BigInteger::Square(vector<int32_t> &w, const vector<int32_t> &x)
{
	uint64_t c, v, prod;
	int32_t wBase, i, j;

	// Note: this method allows w to be only (2 * x.Length - 1) words if result will fit
	// if (w.Length != 2 * x.Length)
	// throw new ArgumentException("no I don't think so...");

	wBase = w.size() - 1;

	i = x.size() - 1;

	while (i > 0)
	{
		v = uint32_t(x[i]);

		c = v * v + uint32_t(w[wBase]);
		w[wBase] = int32_t(c);
		c = c >> 32;

		j = i - 1;
		while (j >= 0)
		{
			prod = v * uint32_t(x[j]);

			wBase--;
			c = c + ((uint32_t(w[wBase]) & UIMASK) + (uint32_t(prod) << 1));
			w[wBase] = int32_t(c);
			c = (c >> 32) + (prod >> 31);
			j--;
		}

		wBase--;
		c = c + uint32_t(w[wBase]);
		w[wBase] = int32_t(c);

		wBase--;
		if (wBase >= 0) w[wBase] = int32_t(c >> 32);

		wBase = wBase + i;
		i--;
	}

	c = uint32_t(x[0]);

	c = (c * c) + uint32_t(w[wBase]);
	w[wBase] = int32_t(c);

	wBase--;
	if (wBase >= 0) w[wBase] = w[wBase] + int32_t(c >> 32);

	return w;
}

vector<int32_t> BigInteger::Multiply(vector<int32_t> &x, const vector<int32_t> &y, const vector<int32_t> &z)
{
	int32_t i, xBase, j;
	int64_t a, val;

	i = z.size();

	if (i < 1) return x;

	xBase = x.size() - y.size();

	do
	{
		i--;
		a = z[i] & IMASK;
		val = 0;

		if (a != 0)
		{
			j = y.size() - 1;

			while (j >= 0)
			{
				val = val + (a * (y[j] & IMASK) + (x[xBase + j] & IMASK));

				x[xBase + j] = int32_t(val);

				val = int64_t(uint64_t(val) >> 32);
				j--;
			}
		}

		xBase--;

		if (xBase >= 0) x[xBase] = int32_t(val);

	} while (i > 0);

	return x;
}

void BigInteger::MontgomeryReduce(vector<int32_t> &x, const vector<int32_t> &m, const uint32_t mDash)
{
	int32_t n, i, j;
	uint32_t x0;
	uint64_t t, carry;

	// NOTE: Not a general purpose reduction (which would allow x up to twice the bitlength of m)
	n = m.size();

	i = n - 1;

	while (i >= 0)
	{
		x0 = uint32_t(x[n - 1]);
		t = uint32_t(x0 * mDash);

		carry = t * uint32_t(m[n - 1]) + x0;

		carry = carry >> 32;

		j = n - 2;

		while (j >= 0)
		{
			carry = carry + (t * uint32_t(m[j]) + uint32_t(x[j]));
			x[j + 1] = int32_t(carry);
			carry = carry >> 32;
			j--;
		}

		x[0] = int32_t(carry);
		i--;
	}

	if (CompareTo(0, x, 0, m) >= 0)
		Subtract(0, x, 0, m);

}

void BigInteger::MultiplyMonty(vector<int32_t> &a, vector<int32_t> &x, const vector<int32_t> &y,
	const vector<int32_t> &m, const uint32_t mDash, const bool smallMontyModulus)
{
	int32_t n, aMax, j, i;
	uint32_t a0, y0;
	uint64_t carry, t, prod1, prod2, xi;

	n = m.size();

	if (n == 1)
	{
		x[0] = int32_t(MultiplyMontyNIsOne(uint32_t(x[0]), uint32_t(y[0]), uint32_t(m[0]), mDash));

		return;
	}

	y0 = uint32_t(y[n - 1]);

	xi = uint32_t(x[n - 1]);

	carry = xi * y0;

	t = uint32_t(uint32_t(carry) * mDash);

	prod2 = t * uint32_t(m[n - 1]);
	carry = carry + uint32_t(prod2);
	carry = (carry >> 32) + (prod2 >> 32);

	j = n - 2;
	while (j >= 0)
	{
		prod1 = xi * uint32_t(y[j]);
		prod2 = t * uint32_t(m[j]);

		carry = carry + ((prod1 & UIMASK) + uint32_t(prod2));
		a[j + 2] = int32_t(carry);
		carry = (carry >> 32) + (prod1 >> 32) + (prod2 >> 32);
		j--;
	}

	a[1] = int32_t(carry);
	aMax = int32_t(carry >> 32);

	i = n - 2;
	while (i >= 0)
	{
		a0 = uint32_t(a[n]);
		xi = uint32_t(x[i]);

		prod1 = xi * y0;
		carry = (prod1 & UIMASK) + a0;
		t = uint32_t(uint32_t(carry) * mDash);

		prod2 = t * uint32_t(m[n - 1]);
		carry = carry + uint32_t(prod2);
		carry = (carry >> 32) + (prod1 >> 32) + (prod2 >> 32);

		j = n - 2;
		while (j >= 0)
		{
			prod1 = xi * uint32_t(y[j]);
			prod2 = t * uint32_t(m[j]);

			carry = carry + ((prod1 & UIMASK) + uint32_t(prod2) + uint32_t(a[j + 1]));
			a[j + 2] = int32_t(carry);
			carry = (carry >> 32) + (prod1 >> 32) + (prod2 >> 32);
			j--;
		}

		carry = carry + uint32_t(aMax);
		a[1] = int32_t(carry);
		aMax = int32_t(carry >> 32);
		i--;
	}

	a[0] = aMax;

	if ((!smallMontyModulus) && (CompareTo(0, a, 0, m) >= 0))
		Subtract(0, a, 0, m);

	memmove(&x[0], &a[1], n * sizeof(int32_t));

}

void BigInteger::SquareMonty(vector<int32_t> &a, vector<int32_t> &x, const vector<int32_t> &m,
	const uint32_t mDash, const bool smallMontyModulus)
{
	int32_t n, aMax, j, i;
	uint32_t xVal, a0;
	uint64_t x0, carry, t, prod1, prod2, xi;

	n = m.size();

	if (n == 1)
	{
		xVal = uint32_t(x[0]);
		x[0] = int32_t(MultiplyMontyNIsOne(xVal, xVal, uint32_t(m[0]), mDash));

		return;
	}

	x0 = uint32_t(x[n - 1]);

	carry = x0 * x0;

	t = uint32_t(uint32_t(carry) * mDash);

	prod2 = t * uint32_t(m[n - 1]);
	carry = carry + uint32_t(prod2);

	carry = (carry >> 32) + (prod2 >> 32);

	j = n - 2;
	while (j >= 0)
	{
		prod1 = x0 * uint32_t(x[j]);
		prod2 = t * uint32_t(m[j]);

		carry = carry + ((prod2 & UIMASK) + (uint32_t(prod1) << 1));
		a[j + 2] = int32_t(carry);
		carry = (carry >> 32) + (prod1 >> 31) + (prod2 >> 32);
		j--;
	}

	a[1] = int32_t(carry);
	aMax = int32_t(carry >> 32);

	i = n - 2;
	while (i >= 0)
	{
		a0 = uint32_t(a[n]);
		t = uint32_t(a0 * mDash);

		carry = t * uint32_t(m[n - 1]) + a0;
		carry = carry >> 32;

		j = n - 2;
		while (j > i)
		{
			carry = carry + (t * uint32_t(m[j]) + uint32_t(a[j + 1]));
			a[j + 2] = int32_t(carry);
			carry = carry >> 32;
			j--;
		}

		xi = uint32_t(x[i]);

		prod1 = xi * xi;
		prod2 = t * uint32_t(m[i]);

		carry = carry + ((prod1 & UIMASK) + uint32_t(prod2) + uint32_t(a[i + 1]));
		a[i + 2] = int32_t(carry);
		carry = (carry >> 32) + (prod1 >> 32) + (prod2 >> 32);

		j = i - 1;
		while (j >= 0)
		{
			prod1 = xi * uint32_t(x[j]);
			prod2 = t * uint32_t(m[j]);

			carry = carry + ((prod2 & UIMASK) + (uint32_t(prod1) << 1) + uint32_t(a[j + 1]));
			a[j + 2] = int32_t(carry);
			carry = (carry >> 32) + (prod1 >> 31) + (prod2 >> 32);
			j--;
		}

		carry = carry + uint32_t(aMax);
		a[1] = int32_t(carry);
		aMax = int32_t(carry >> 32);
		i--;
	}

	a[0] = aMax;

	if ((!smallMontyModulus) && (CompareTo(0, a, 0, m) >= 0))
		Subtract(0, a, 0, m);

	memmove(&x[0], &a[1], n * sizeof(int32_t));

}

uint32_t BigInteger::MultiplyMontyNIsOne(const uint32_t x, const uint32_t y, const uint32_t m, const uint32_t mDash)
{
	uint64_t carry, um, prod2;
	uint32_t t;

	carry = uint64_t(uint64_t(x) * y);
	t = uint32_t(uint32_t(carry) * mDash);
	um = m;
	prod2 = um * t;
	carry = carry + uint32_t(prod2);

	carry = (carry >> 32) + (prod2 >> 32);
	if (carry > um) carry = carry - um;

	return uint32_t(carry);
}

vector<int32_t> BigInteger::ShiftLeft(const vector<int32_t> &mag, const int32_t n)
{
	int32_t nInts, nBits, magLen, i, nBits2, highBits, m, j, Next;
	vector<int32_t> newMag;

	nInts = int32_t(uint32_t(n) >> 5);
	nBits = n & 0x1F;
	magLen = mag.size();

	if (nBits == 0)
	{
		newMag = vector<int32_t>(magLen + nInts);
		memmove(&newMag[0], &mag[0], mag.size() * sizeof(int32_t));
	}
	else
	{
		i = 0;
		nBits2 = 32 - nBits;
		highBits = int32_t(uint32_t(uint32_t(mag[0]) >> nBits2));

		if (highBits != 0)
		{
			newMag = vector<int32_t>(magLen + nInts + 1);

			newMag[i] = highBits;
			i++;
		}
		else
			newMag = vector<int32_t>(magLen + nInts);

		m = mag[0];

		j = 0;
		while (j < (magLen - 1))
		{
			Next = mag[j + 1];

			newMag[i] = (uint32_t(m) << nBits) | int32_t(uint32_t(Next) >> nBits2);
			i++;
			m = Next;
			j++;
		}

		newMag[i] = int32_t(uint32_t(mag[magLen - 1]) << nBits);
	}

	return newMag;
}

void BigInteger::ShiftRightInPlace(const int32_t start, vector<int32_t> &mag, const int32_t n)
{
	int32_t nInts, nBits, magEnd, delta, i, nBits2, m, Next;

	nInts = int32_t(uint32_t(n) >> 5) + start;
	nBits = n & 0x1F;
	magEnd = mag.size() - 1;

	if (nInts != start)
	{
		delta = (nInts - start);

		i = magEnd;
		while (i >= nInts)
		{
			mag[i] = mag[i - delta];
			i--;
		}

		i = nInts - 1;
		while (i >= start)
		{
			mag[i] = 0;
			i--;
		}

	}

	if (nBits != 0)
	{
		nBits2 = 32 - nBits;
		m = mag[magEnd];

		i = magEnd;
		while (i > nInts)
		{
			Next = mag[i - 1];

			mag[i] = int32_t(uint32_t(m) >> nBits) | (Next << nBits2);
			m = Next;
			i--;
		}

		mag[nInts] = int32_t(uint32_t(mag[nInts]) >> nBits);
	}

}

void BigInteger::ShiftRightOneInPlace(const int32_t start, vector<int32_t> &mag)
{
	int32_t i, m, Next;

	i = mag.size();
	m = mag[i - 1];

	i--;
	while (i > start)
	{
		Next = mag[i - 1];
		mag[i] = (int32_t(uint32_t(m) >> 1)) | (Next << 31);
		m = Next;
		i--;
	}

	mag[start] = int32_t(uint32_t(mag[start]) >> 1);
}

vector<int32_t> BigInteger::Subtract(const int32_t xStart, vector<int32_t> &x, const int32_t yStart, const vector<int32_t> &y)
{
	int32_t iT, iV, borrow;
	int64_t m;

	iT = x.size();
	iV = y.size();

	borrow = 0;

	do
	{
		iT--;
		iV--;
		m = (uint64_t(x[iT]) & IMASK) - ((uint64_t(y[iV]) & IMASK) + borrow);
		// fixed precedence bug :)
		x[iT] = int32_t(m);
		borrow = int32_t(uint64_t(m) >> 63);
	} while (iV > yStart);

	if (borrow != 0)
	{
		iT--;
		x[iT] = x[iT] - 1;
		while (x[iT] == -1)
		{
			iT--;
			x[iT] = x[iT] - 1;
		}
	}

	return x;
}

vector<int32_t> BigInteger::doSubBigLil(const vector<int32_t> &bigMag, const vector<int32_t> &lilMag)
{
	vector<int32_t> res;

	res = bigMag;

	return Subtract(0, res, 0, lilMag);
}

void BigInteger::AppendZeroExtendedString(ostringstream &sl, const string s, const int32_t minLength)
{
	int32_t len;

	len = s.size();
	while (len < minLength)
	{
		sl << '0';
		len++;
	}

	sl << s;
}

BigInteger BigInteger::CreateUValueOf(const uint64_t value)
{
	int32_t msw, lsw;

	msw = int32_t(value >> 32);
	lsw = int32_t(value);

	if (msw != 0)
		return BigInteger(1, vector<int32_t>({ msw, lsw }), false);

	if (lsw != 0)
	{
		BigInteger result = BigInteger(1, vector<int32_t>({ lsw }), false);
		// Check for a power of two

		if ((lsw & -lsw) == lsw) result.nBits = 1;

		return result;
	}

	return Zero;
}

BigInteger BigInteger::CreateValueOf(const int64_t value)
{
	if (value < 0)
	{
		if (value == INT64_MIN)
			return CreateValueOf(~value).Not();

		return CreateValueOf(-value).Negate();
	}

	return CreateUValueOf(uint64_t(value));
}

string BigInteger::IntToBin(const int32_t _input)
{
	int32_t input = _input;
	vector<char> bits;

	while (input != 0)
	{
		if ((input & 1) == 1) bits.push_back('1');
		else
			bits.push_back('0');

		input = uint32_t(input) >> 1;
	}

	string result(bits.begin(), bits.end());

	reverse(result.begin(), result.end());

	return result;
}

string BigInteger::IntToOctal(const int32_t _input)
{
	int32_t input = _input;
	vector<char> bits;
	char temp[25];

	while (input != 0)
	{
		ui64toa_s(uint64_t(input & 7), &(temp[0]), 25, 10);

		bits.push_back(temp[0]);
		input = uint32_t(input) >> 3;
	}

	string result(bits.begin(), bits.end());

	reverse(result.begin(), result.end());

	return result;
}

class Point {
public: 
	BigInteger x, y;
	Point(BigInteger x, BigInteger y) : x(x), y(y) {}
	bool operator == (const Point & other) const {
		return x.Equals(other.x) && y.Equals(other.y);
	}
	bool operator != (const Point & other) const {
		return !(*this == other);
	}
};

class EllipticCurve {
private: 
	BigInteger a, b, p;
	BigInteger bigZero = BigInteger("0");
	BigInteger big3 = BigInteger("3");
	BigInteger big2 = BigInteger("2");
	BigInteger big1 = BigInteger("1");

public:
	EllipticCurve(BigInteger a, BigInteger b, BigInteger p) : a(a), b(b), p(p)
	{
		BigInteger big4("4");
		BigInteger big27("27");

		//4 * a * a * a
		//27 * b * b
		BigInteger big4aaa = big4.Multiply(a).Multiply(a).Multiply(a);
		BigInteger big27bb = big27.Multiply(b).Multiply(b);
		if (big4aaa.Subtract(big27bb).Mod(p).Equals(bigZero))
		{
			throw std::invalid_argument("The given parameters do not form an elliptic curve.");
		}
	}

	Point pointAddition(const Point & P, const Point & Q) const
	{
		if (P == Point(bigZero, bigZero)) return Q;
		if (Q == Point(bigZero, bigZero)) return P;
		if (P.x.Equals(Q.x) && !P.y.Equals(Q.y)) return Point(bigZero, bigZero);
		BigInteger m;

		if (P == Q)
		{
			//m = (3 * P.x * P.x + a) * modInverse(2 * P.y, p) % p;

			BigInteger big3PA = big3.Multiply(P.x).Multiply(P.x).Add(a);
			BigInteger bigModInvere = big2.Multiply(P.y).ModInverse(p);
			m = big3PA.Multiply(bigModInvere).Mod(p);
		}
		else
		{
			//m = (Q.y - P.y) * modInverse(Q.x - P.x, p) % p;
			BigInteger bigQySPy = Q.y.Subtract(P.y);
			BigInteger bigModInvere = Q.x.Subtract(P.x).ModInverse(p);
			m = bigQySPy.Multiply(bigModInvere).Mod(p);
		}

		BigInteger x_r = (m.Square().Subtract(P.x).Subtract(Q.x)).Mod(p);
		BigInteger y_r = (m.Multiply(P.x.Subtract(x_r)).Subtract(P.y)).Mod(p);
		if (x_r.CompareTo(bigZero) < 0) x_r = x_r.Add(p);
		if (y_r.CompareTo(bigZero) < 0) y_r = y_r.Add(p);
		return Point(x_r, y_r);
	}
};

std::string toUpperCase(const std::string& input)
{
	std::string result = input;
	std::transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}

int main(int argc, const char * argv[])
{
	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " test.inp test.out" << std::endl;
		return 1;
	}

	try
	{
		std::string inputFileName = std::string(argv[1]);
		std::string outputFileName = std::string(argv[2]);
		std::ifstream inputFile(inputFileName);
		std::ofstream outputFile(outputFileName);
		if (!inputFile.is_open() || !outputFile.is_open())
		{
			std::cerr << "Error opening files!" << std::endl;
			return 1;
		}

		std::string p_hex, a_hex, b_hex, xP_hex, yP_hex, xQ_hex, yQ_hex;
		std::getline(inputFile, p_hex, ' ');
		std::getline(inputFile, a_hex, ' ');
		std::getline(inputFile, b_hex);
		std::getline(inputFile, xP_hex, ' ');
		std::getline(inputFile, yP_hex);
		std::getline(inputFile, xQ_hex, ' ');
		std::getline(inputFile, yQ_hex);

		std::cout << "p: " << p_hex << std::endl;
		std::cout << "a: " << a_hex << std::endl;
		std::cout << "b: " << b_hex << std::endl;
		std::cout << "xP: " << xP_hex << std::endl;
		std::cout << "yP: " << yP_hex << std::endl;
		std::cout << "xQ: " << xQ_hex << std::endl;
		std::cout << "yQ: " << yQ_hex << std::endl;

		BigInteger::Boot();

		BigInteger p(p_hex, 16);
		BigInteger a(a_hex, 16);
		BigInteger b(b_hex, 16);
		BigInteger xP(xP_hex, 16);
		BigInteger yP(yP_hex, 16);
		BigInteger xQ(xQ_hex, 16);
		BigInteger yQ(yQ_hex, 16);

		EllipticCurve curve(a, b, p);

		Point P(xP, yP);
		Point Q(xQ, yQ);

		Point R = curve.pointAddition(P, Q);

		std::cout << "Result. " << std::endl;
		std::cout << "xR: " << toUpperCase(R.x.ToString(16)) << std::endl;
		std::cout << "yR: " << toUpperCase(R.y.ToString(16)) << std::endl;

		outputFile << toUpperCase(R.x.ToString(16)) << ' ' << toUpperCase(R.y.ToString(16));

		inputFile.close();
		outputFile.close();

		std::cout << "Done writing to file." << std::endl;
		return 0;
	}
	catch (std::exception& ex)
	{
		std::cerr << "Exception: " << ex.what() << std::endl;
		return 0;
	}
}

#endif
