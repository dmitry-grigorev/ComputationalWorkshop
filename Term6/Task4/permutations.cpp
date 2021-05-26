#include "pch.h"
#include <iostream>
#include <random>
#include "permutations.h"

void Permutation::generate()
{
	int j = 0, buf = 0;
	for (int i = permlen - 1; i > 0; --i)
	{
		j = gen->randint(0, i);

		buf = permutation[i];
		permutation[i] = permutation[j];
		permutation[j] = buf;
	}

	CaclCyclesLensDist();
}

void Permutation::CaclCyclesLensDist()
{
	int* buf = new int[permlen];
	int currlen = 0;

	for (int i = 0; i < permlen; ++i)
	{
		buf[i] = permutation[i];
		cycleslensdist[i] = 0;
	}

	for (int i = 0; i < permlen; ++i)
	{
		currlen = 0;
		int k = i, kcopy = k;
		while (buf[k] != 0)
		{
			++currlen;
			kcopy = k;
			k = buf[k] - 1;
			buf[kcopy] = 0;
		}

		if (currlen > 0)
			++cycleslensdist[currlen - 1];
	}

	delete buf;
}