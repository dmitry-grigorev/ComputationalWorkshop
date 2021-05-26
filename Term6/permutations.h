#pragma once
class Generator
{
	unsigned int iu;
	const double flt = 0.232830643654e-9;
protected:
	std::mt19937 *gen;
public:
	Generator(const unsigned long seed = 0) {
		gen = new std::mt19937(seed);
	};
private:
	double rnunif(void)
	{
		return (flt* (*gen)());
	}
public:
	int randint(const int low, const int high)
	{
		return low + int((high + 1)*rnunif());
	}

	~Generator()
	{
		delete gen;
	}
};

class Permutation
{
	Generator *gen;
	int* permutation;
	unsigned int permlen;
	int* cycleslensdist;
public:
	Permutation(const unsigned int permlen = 10, const unsigned long seed = 0)
	{
		gen = new Generator(seed);
		permutation = new int[permlen];
		this->permlen = permlen;
		for (int i = 0; i < permlen; ++i)
			permutation[i] = i + 1;

		cycleslensdist = new int[permlen] {};
	}

	~Permutation()
	{
		delete gen;
		delete permutation;
		delete cycleslensdist;
	}

	int operator [](const int i)
	{
		return permutation[i];
	}

	void printPermutation()
	{
		for (int i = 0; i < permlen; ++i)
		{
			std::cout << permutation[i] << ' ';
		}
		std::cout << std::endl;
	}

	void printCyclesLensDist()
	{
		for (int i = 0; i < permlen; ++i)
		{
			std::cout << cycleslensdist[i] << ' ';
		}
		std::cout << std::endl;
	}

	void generate();

	void CaclCyclesLensDist();

	int* getCyclesLensDist()
	{
		return cycleslensdist;
	}
};
