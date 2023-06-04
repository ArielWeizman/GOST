#include <string> 
#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <ctime>

using namespace std;

#define NGOODPAIRS 256
#define BLOCKSIZE 256
#define MAXMESSAGELEN 64

int DS5 = 0;
int DS4 = 0;
int DS1 = 0;
int DS2 = 0;

uint64_t ivalue;
int len = 1;

short s_block[8][16] = {};
//Bank S-boxes
short bankS[8][16] = { { 4, 0xA, 9, 2, 0xD, 8, 0, 0xE, 6, 0xB, 1, 0xC, 7, 0xF, 5, 3 },
{ 0xE, 0xB, 4, 0xC, 6, 0xD, 0xF, 0xA, 2, 3, 8, 1, 0, 7, 5, 9 },
{ 5, 8, 1, 0xD, 0xA, 3, 4, 2, 0xE, 0xF, 0xC, 7, 6, 0, 9, 0xB },
{ 7, 0xD, 0xA, 1, 0, 8, 9, 0xF, 0xE, 4, 6, 0xC, 0xB, 2, 5, 3 },
{ 6, 0xC, 7, 1, 5, 0xF, 0xD, 8, 4, 0xA, 9, 0xE, 0, 3, 0xB, 2 },
{ 4, 0xB, 0xA, 0, 7, 2, 1, 0xD, 3, 6, 8, 5, 9, 0xC, 0xF, 0xE },
{ 0xD, 0xB, 4, 1, 3, 0xF, 5, 9, 0, 0xA, 0xE, 7, 6, 8, 2, 0xC },
{ 1, 0xF, 0xD, 0, 5, 7, 0xA, 4, 9, 2, 3, 0xE, 6, 0xB, 8, 0xC } };

short sInverse[8][16];

//gost2 sboxes
short gost2S[8][16] = { { 6,10,15,4,3,8,5,0,13,14,7,1,2,11,12,9 },
{ 6,10,15,4,3,8,5,0,13,14,7,1,2,11,12,9 },
{ 6,10,15,4,3,8,5,0,13,14,7,1,2,11,12,9 },
{ 6,10,15,4,3,8,5,0,13,14,7,1,2,11,12,9 },
{ 14,0,8,1,7,10,5,6,13,2,4,9,3,15,12,11 },
{ 14,0,8,1,7,10,5,6,13,2,4,9,3,15,12,11 },
{ 14,0,8,1,7,10,5,6,13,2,4,9,3,15,12,11 },
{ 14,0,8,1,7,10,5,6,13,2,4,9,3,15,12,11 } };

void initSInv()
{
	for (int s = 0; s < 8; s++)
		for (int v = 0; v < 16; v++)
			sInverse[s][s_block[s][v]] = v;
}

short ddt0[16][16] = {};
short ddt2[16][16] = {};
short ddt3[16][16] = {};
short ddt5[16][16] = {};

void initDDT()
{
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt0[i ^ j][s_block[0][i] ^ s_block[0][j]]++;
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt2[i ^ j][s_block[2][i] ^ s_block[2][j]]++;
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt3[i ^ j][s_block[3][i] ^ s_block[3][j]]++;
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt5[i ^ j][s_block[5][i] ^ s_block[5][j]]++;
}

void calcDDT(short sb[16], int ddt[16][16])
{
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt[i ^ j][sb[i] ^ sb[j]]++;
}

void printDDT(short sbox[16])
{
	ofstream file;
	file.open("ddt.txt");
	int ddt[16][16] = {};
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt[i ^ j][sbox[i] ^ sbox[j]]++;
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			file << ddt[i][j] << "\t";
		}
		file << endl;
	}
	file.close();
}

void sort2dimArrayByACol(int A[][4], int n, int col)
{
	for (int i = 1; i < n; i++)
	{
		for (int j = i - 1; j >= 0 && A[j + 1][col] < A[j][col]; j--)
		{
			//Swap the values in all the  rows
			for (int r = 0; r < 4; r++)
			{
				swap(A[j][r], A[j + 1][r]);
			}
		}
	}
}

void Afunction(const uint8_t M[32], uint8_t res[32])
{
	uint64_t temp1[4] = {};
	uint64_t temp2[4] = {};
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 8; j++)
			temp1[i] |= ((uint64_t)(M[8 * i + j]) << 8 * j);

	for (int i = 0; i < 3; i++)
		temp2[i] = temp1[i + 1];
	temp2[3] = (temp1[0] ^ temp1[1]);
	for (int i = 0; i < 32; i++)
		res[i] = (uint8_t)(temp2[i / 8] >> (8 * (i % 8)));

}

int phiFunc(int val)
{
	int i = (val % 4);
	int k = (val / 4);
	int res = (8 * i + k);
	return res;
}

void PFunc(uint8_t val[32])
{
	uint8_t temp[32] = {};
	for (int i = 0; i < 32; i++)
		temp[i] = val[phiFunc(i)];
	for (int i = 0; i < 32; i++)
		val[i] = temp[i];
}

void randomShuffle(int arr[16], int arrSize)
{
	for (int i = 0; i < arrSize; i++)
	{
		int index = rand() % (arrSize - i);
		swap(arr[index], arr[arrSize - i - 1]);
	}
}

uint64_t round(uint64_t _block, uint32_t subkey, short s_block[][16])
{
	uint64_t block = 0;
	uint32_t right = _block;
	uint32_t left;
	uint32_t N;
	uint32_t SN = 0;
	uint32_t right1 = right;
	left = _block >> 32;
	N = (right + subkey) % 4294967296;


	for (int j = 0; j <= 7; j++)
	{
		uint8_t Ni = (N >> (4 * j)) % 16;
		Ni = s_block[j][Ni]; // substitution through optionInd-blocks.

		uint32_t mask = Ni;
		mask <<= (4 * j);

		SN = SN | mask;
	}
	N = SN;

	uint32_t mask = N << 11;
	N = (N >> 21) | mask;

	right = N ^ left;
	left = right1;
	block = block | left;
	block = block << 32;
	block = block | right;

	return block;
}

uint32_t NEcr = 0;
uint64_t encrypt(uint64_t _block, uint32_t *key, short s_block[][16], int nRounds)
{
	NEcr++;
	uint64_t block = _block;
	uint32_t left = 0;

	for (int rnd = 0; rnd < nRounds; rnd++)
	{
		if (rnd < 24)
			block = round(block, key[rnd % 8], s_block);
		else
			block = round(block, key[31 - rnd], s_block);
	}

	block = (block << 32) | (block >> 32);

	return block;
}

uint64_t encrypt2(uint64_t _block, uint32_t *key, short s_block[][16], int nRounds)
{
	uint64_t block = _block;
	uint32_t left = 0;

	for (int rnd = 0; rnd < nRounds; rnd++)
	{
		if (rnd < 8)
			block = round(block, key[rnd % 8], s_block);
		else if (rnd < 16)
			block = round(block, key[(rnd + 3) % 8], s_block);
		else if (rnd < 24)
			block = round(block, key[(rnd + 5) % 8], s_block);
		else
			block = round(block, key[(38 - rnd) % 8], s_block);
	}

	return block;
}

uint64_t decrypt(uint64_t _block, uint32_t *key, short s_block[][16])
{
	uint64_t block = _block;
	uint32_t left = 0;

	for (int i = 0; i <= 7; i++)
		block = round(block, key[i], s_block);

	for (int k = 1; k <= 3; k++)
		for (int i = 7; i >= 0; i--)
			block = round(block, key[i], s_block);

	block = (block << 32) | (block >> 32);


	return block;
}

void buildInValOutDiffs(int s, uint64_t c1s[], uint64_t c2s[],
	uint32_t ke, int inValsOutDiffs[][4], bool moreC, uint32_t rk, int relatedBit)
{
	for (int i = 0; i < NGOODPAIRS + moreC * NGOODPAIRS; i++)
	{
		uint32_t c1PlusK = c1s[i] + ke;
		uint32_t c2PlusK = c2s[i] + (ke ^ ((uint32_t)rk << (7 + relatedBit)));
		inValsOutDiffs[i][0] = ((c1PlusK & (uint32_t)(0xf << 4 * s)) >> 4 * s);
		inValsOutDiffs[i][1] = ((c2PlusK & (uint32_t)(0xf << 4 * s)) >> 4 * s);
		inValsOutDiffs[i][3] = i;
		for (int b = 0; b < 4; b++)
		{
			int pos = 32 + ((27 + b + 4 * (s - 4)) % 32);
			uint64_t mask = ((uint64_t)1 << pos);
			uint64_t c1Bit = (c1s[i] & mask);
			uint64_t c2Bit = (c2s[i] & mask);
			inValsOutDiffs[i][2] |= (((c1s[i] & mask) ^ (c2s[i] & mask))
				>> (pos - b));
		}
	}

	sort2dimArrayByACol(inValsOutDiffs, NGOODPAIRS + moreC * NGOODPAIRS, 0);
}

bool isWrongKey(int inValsOutDiffs[][4], int inds[17], bool moreC)
{
	int ind = 0;
	while (inValsOutDiffs[inds[ind]][0] == ind && ind<15)
	{
		inds[ind + 1]++;
		if (inValsOutDiffs[inds[ind + 1]][0] == ind + 1)
		{
			ind++;
			inds[ind + 1] = inds[ind];
		}
		else if (inValsOutDiffs[inds[ind + 1]][0] > ind + 1)
		{
			ind++;
			inds[ind + 1] = inds[ind];
		}
	}
	inds[16] = NGOODPAIRS + moreC * NGOODPAIRS;

	//if (probl)
	//{
	//	for (int t = 0; t < 17; t++)
	//		cout << inds[t] << endl;
	//}
	//Check consistency
	//for (int i1 = 0; i1 < NGOODPAIRS + moreC * NGOODPAIRS - 1; i1++)
	//{
	//	for (int i2 = i1 + 1; i2 < NGOODPAIRS + moreC * NGOODPAIRS; i2++)
	//	{
	//		if (inValsOutDiffs[i1][0] == inValsOutDiffs[i2][0])
	//		{
	//			if ((inValsOutDiffs[i1][1] == inValsOutDiffs[i2][1])
	//				&& (inValsOutDiffs[i1][2] != inValsOutDiffs[i2][2]))
	//				return true;
	//		}
	//	}
	//}
	for (int ii = 0; ii < 16; ii++)
	{
		for (int i1 = inds[ii]; i1 < inds[ii + 1] - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < inds[ii + 1]; i2++)
			{
				if ((inValsOutDiffs[i1][1] == inValsOutDiffs[i2][1]) &&
					(inValsOutDiffs[i1][2] != inValsOutDiffs[i2][2]))
					return true;
			}
		}
	}
	return false;
}

void buildSbox(int rSbox[8][16], int inValsOutDiffs[][4], int inds[17], int s, bool moreC)
{
	//for (int v = 0; v < inds[1]; v++)
	//	rSbox[s][inValsOutDiffs[v][1]] = inValsOutDiffs[v][2];
	vector<bool> done;
	done.resize(NGOODPAIRS + moreC * NGOODPAIRS);
	for (int i = 0; i < 15; i++)
	{
		bool endB = true;
		for (int j = 1; j < 16; j++)
		{
			if (!rSbox[s][j])
			{
				endB = false;
				break;
			}
		}
		if (endB)
			break;
		for (int v = 0; v < NGOODPAIRS + moreC * NGOODPAIRS; v++)
		{
			if (done[v])
				continue;
			if (!inValsOutDiffs[v][0] && inValsOutDiffs[v][1])
			{
				rSbox[s][inValsOutDiffs[v][1]] = inValsOutDiffs[v][2];
				done[v] = true;
			}
			else if (!inValsOutDiffs[v][1] && inValsOutDiffs[v][0])
			{
				rSbox[s][inValsOutDiffs[v][0]] = inValsOutDiffs[v][2];
				done[v] = true;
			}
			else if (rSbox[s][inValsOutDiffs[v][1]])
			{
				rSbox[s][inValsOutDiffs[v][0]] = rSbox[s][inValsOutDiffs[v][1]]
					^ inValsOutDiffs[v][2];
				done[v] = true;
			}
			else if (rSbox[s][inValsOutDiffs[v][0]])
			{
				rSbox[s][inValsOutDiffs[v][1]] = rSbox[s][inValsOutDiffs[v][0]]
					^ inValsOutDiffs[v][2];
				done[v] = true;
			}
		}
	}

	//for (int v = 0; v < inds[1]; v++)
	//	rSbox[s][inValsOutDiffs[v][1]] = inValsOutDiffs[v][2];

	//for (int i = 0; i < 15; i++)
	//{
	//	for (int sVal = 1; sVal < 16; sVal++)
	//	{
	//		if (rSbox[s][sVal])
	//			continue;
	//		for (int inSVal = inds[sVal]; inSVal < inds[sVal + 1]; inSVal++)
	//		{
	//			if (rSbox[s][inValsOutDiffs[inSVal][1]])
	//			{
	//				rSbox[s][sVal] = rSbox[s][inValsOutDiffs[inSVal][1]]
	//					^ inValsOutDiffs[inSVal][2];
	//				break;
	//			}
	//		}
	//	}
	//}

	int remains = 0;
	bool empties[16] = {};
	empties[0] = true;
	bool Res[16] = {};
	Res[0] = true;
	int empt;
	for (int sVal = 1; sVal < 16; sVal++)
	{
		if (!rSbox[s][sVal])
		{
			empt = sVal;
			remains++;
		}
		else
		{
			empties[sVal] = true;
			Res[rSbox[s][sVal]] = true;
		}
	}
	if (remains == 1)
	{
		for (int sVal = 1; sVal < 16; sVal++)
		{
			if (!Res[sVal])
				rSbox[s][empt] = sVal;
		}
	}
}

bool s1isBelong(int s1Options[16][16], int s2[8][16])
{
	for (int optionInd = 0; optionInd < 16; optionInd++)
	{
		bool isBelong = true;
		for (int v = 0; v < 16; v++)
			if (s1Options[optionInd][v] != s2[1][v])
				isBelong = false;
		if (isBelong)
			return true;
	}
	return false;
}

int SUC = 0;
int SUC1 = 0;
int SUC2 = 0;
int SUCSUC = 0;
int SUCSUCOI = 0;

void checkfrst16bitK0(uint64_t c1s[NGOODPAIRS], uint64_t c2s[NGOODPAIRS], bool goodKeys[65536], uint32_t theKey)
{
	bool relevantC[NGOODPAIRS] = {};
	int NRelevantC = 0;
	for (int i = 0; i < NGOODPAIRS; i++)
	{
		if ((c1s[i] & 0x8000) ^ (c2s[i] & 0x8000))
		{
			relevantC[i] = true;
			NRelevantC++;
		}
	}
	vector<uint64_t> rc1s(NRelevantC);
	vector<uint64_t> rc2s(NRelevantC);
	int ii = 0;
	for (int i = 0; i < NGOODPAIRS; i++)
	{
		if (relevantC[i])
		{
			rc1s[ii] = c1s[i];
			rc2s[ii] = c2s[i];
			ii++;
		}
	}

	for (int k = 0; k < 65536; k++)
	{
		vector<int[4]> inValsOutDiffs(NRelevantC);

		if (k == (theKey & 0xffff))
			int gdsgsd = 0;

		for (int i = 0; i < NRelevantC; i++)
		{
			uint32_t c1PlusK = rc1s[i] + k;
			uint32_t c2PlusK = rc2s[i] + k;
			inValsOutDiffs[i][0] = ((c1PlusK & (uint32_t)(0xf << 12)) >> 12);
			inValsOutDiffs[i][1] = ((c2PlusK & (uint32_t)(0xf << 12)) >> 12);
			inValsOutDiffs[i][3] = i;
			for (int b = 0; b < 4; b++)
			{
				int pos = 55 + b;
				uint64_t mask = ((uint64_t)1 << pos);
				uint64_t c1Bit = (rc1s[i] & mask);
				uint64_t c2Bit = (rc2s[i] & mask);
				inValsOutDiffs[i][2] |= (((rc1s[i] & mask) ^ (rc2s[i] & mask))
					>> (pos - b));
			}
		}

		//for (int i = 1; i < NRelevantC; i++)
		//{
		//	for (int j = i - 1; j >= 0 && inValsOutDiffs[j + 1][0] < inValsOutDiffs[j][0]; j--)
		//	{
		//		//Swap the values in all the  rows
		//		for (int r = 0; r < 4; r++)
		//		{
		//			swap(inValsOutDiffs[j][r], inValsOutDiffs[j + 1][r]);
		//		}
		//	}
		//}

		//Check consistency
		for (int i1 = 0; ((i1 < NRelevantC - 1) && (goodKeys[k])); i1++)
		{
			for (int i2 = i1 + 1; i2 < NRelevantC; i2++)
			{
				if ((inValsOutDiffs[i1][1] == inValsOutDiffs[i2][1]) &&
					(inValsOutDiffs[i1][2] != inValsOutDiffs[i2][2]))
				{
					goodKeys[k] = false;
					break;
				}
			}
		}
	}
	//int remK = 0;
	//for (int expInd = 0; expInd < 65536; expInd++)
	//	remK += goodKeys[expInd];
	//return remK;
}

bool checkfrst7bitK1(vector <uint64_t> rc1s, vector <uint64_t> rc2s, bool goodKeys[256],
	short sbox[8][16], uint32_t ke, int NRelevantC)
{
	for (int i = 0; i < 256; i++)
		goodKeys[i] = 1;
	for (int k = 0; k < 256; k++)
	{
		vector<int[3]> inValsOutDiffs(NRelevantC);

		for (int i = 0; i < NRelevantC; i++)
		{
			uint32_t c1PlusK = rc1s[i] + ke;
			uint32_t c2PlusK = rc2s[i] + ke;
			int aftrSbox1 = ((sbox[2][(c1PlusK & 0xf00000) >> 20]) >> 1);
			aftrSbox1 |= ((sbox[1][(c1PlusK & 0xf000000) >> 24]) << 3);
			aftrSbox1 ^= ((rc1s[i] & 0x7f00000000) >> 32);
			int aftrSbox2 = ((sbox[2][(c2PlusK & 0xf00000) >> 20]) >> 1);
			aftrSbox2 |= ((sbox[1][(c2PlusK & 0xf000000) >> 24]) << 3);
			aftrSbox2 ^= ((rc2s[i] & 0x7f00000000) >> 32);
			inValsOutDiffs[i][0] = (((aftrSbox1 + k) & 0x70) >> 4);
			inValsOutDiffs[i][1] = (((aftrSbox2 + k) & 0x70) >> 4) + 8;
			//inValsOutDiffs[i][3] = i;
			uint64_t mask = 0x78000;
			uint64_t c1Val = (rc1s[i] & mask);
			uint64_t c2Val = (rc2s[i] & mask);
			inValsOutDiffs[i][2] = ((c1Val^c2Val) >> 15);
		}

		//for (int i = 1; i < NRelevantC; i++)
		//{
		//	for (int j = i - 1; j >= 0 && inValsOutDiffs[j + 1][0] < inValsOutDiffs[j][0]; j--)
		//	{
		//		//Swap the values in all the  rows
		//		for (int r = 0; r < 4; r++)
		//		{
		//			swap(inValsOutDiffs[j][r], inValsOutDiffs[j + 1][r]);
		//		}
		//	}
		//}

		//Check consistency
		for (int i1 = 0; ((i1 < NRelevantC - 1) && (goodKeys[k])); i1++)
		{
			for (int i2 = i1 + 1; i2 < NRelevantC; i2++)
			{
				if ((inValsOutDiffs[i1][1] == inValsOutDiffs[i2][1]) &&
					(inValsOutDiffs[i1][2] != inValsOutDiffs[i2][2]))
				{
					goodKeys[k] = false;
					break;
				}
			}
		}
	}

	for (int i = 0; i < 128; i++)
		if (goodKeys[i])
			return true;
	return false;
}

bool earlRej2(uint64_t c1s[NGOODPAIRS], uint64_t c2s[NGOODPAIRS], uint32_t k)
{
	//uint32_t ke = (goodKeys[k / 16]) | ((k % 16) << 16);
	int inValsOutDiffs[NGOODPAIRS][4] = {};

	//if (ke == (theKey & 0xffff))
	//	int gdsgsd = 0;
	buildInValOutDiffs(4, c1s, c2s, k, inValsOutDiffs, 0, 0, 0);

	//Check consistency
	for (int i1 = 0; i1 < NGOODPAIRS - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < NGOODPAIRS; i2++)
		{
			if ((inValsOutDiffs[i1][1] == inValsOutDiffs[i2][1]) &&
				(inValsOutDiffs[i1][2] != inValsOutDiffs[i2][2]))
			{
				return false;
			}
		}
	}

	return true;
}

bool SboxRecoveryAttack0(uint64_t c1s[NGOODPAIRS], uint64_t c2s[NGOODPAIRS], uint32_t theKey,
	vector<uint32_t> &resultKeys, vector<uint64_t> &s4, vector<uint64_t> &s5)
{
	uint64_t hardMask = 0x7ff800ff807fff;
	uint64_t harderMask = 0x7ff800fff07fff;
	uint64_t weakMask = 0x7ff80080007fff;
	int sb[16][16] = {};
	int sbInd = 0;
	bool goodKeys[65536];
	for (int i = 0; i < 65536; i++)
		goodKeys[i] = true;
	checkfrst16bitK0(c1s, c2s, goodKeys, theKey);
	for (uint32_t ke = 0; ke < 65536; ke++)
	{
		if (!goodKeys[ke])
			continue;
		for (int cind = 0; cind < NGOODPAIRS; cind++)
		{
			if ((c1s[cind] & (uint64_t)0x7ff800000f7fff) == (c2s[cind] & (uint64_t)0x7ff800000f7fff))
			{
				if ((c1s[cind] & ((uint64_t)(0x7800000000000000))) ^ (c2s[cind] & ((uint64_t)(0x7800000000000000))))
				{
					uint32_t val1 = (ke % 0x10000) + (c1s[cind] & ((uint32_t)0xffff));
					uint32_t val2 = (ke % 0x10000) + (c2s[cind] & ((uint32_t)0xffff));
					if ((val1 >= 65536) && (val2 >= 65536))
					{
						goodKeys[ke] = 0;
						break;
					}
					if ((val1 < 65536) && (val2 < 65536))
					{
						goodKeys[ke] = 0;
						break;
					}
				}
			}
		}
	}

	int ttt = 0;
	for (int i = 0; i < 65536; i++)
		ttt += goodKeys[i];
	vector<uint32_t> optional16BitKeys;
	optional16BitKeys.resize(ttt);
	int ind = 0;
	for (uint32_t i = 0; i < 65536; i++)
		if (goodKeys[i])
			optional16BitKeys[ind++] = i;

	vector<bool> suc20BitKeys;
	suc20BitKeys.resize(ttt * 16);
	for (int i = 0; i < ttt * 16; i++)
		suc20BitKeys[i] = true;
	for (int i = 0; i < ttt; i++)
	{
		for (uint32_t j = 0; j < 16; j++)
		{
			uint32_t ke = optional16BitKeys[i] | (j << 16);
			for (int cind = 0; cind < NGOODPAIRS; cind++)
			{
				if ((c1s[cind] & harderMask) == (c2s[cind] & harderMask))
				{
					if ((c1s[cind] & ((uint64_t)(0x8000000700000000))) ^ (c2s[cind] & ((uint64_t)(0x8000000700000000))))
					{
						uint32_t val1 = (ke % 0x100000) + (c1s[cind] & ((uint32_t)0xfffff));
						uint32_t val2 = (ke % 0x100000) + (c2s[cind] & ((uint32_t)0xfffff));
						if ((val1 >= 1048576) && (val2 >= 1048576))
						{
							suc20BitKeys[16 * i + j] = 0;
							break;
						}
						if ((val1 < 1048576) && (val2 < 1048576))
						{
							suc20BitKeys[16 * i + j] = 0;
							break;
						}
					}
				}
			}
		}
	}

	ttt = 0;
	for (int i = 0; i < optional16BitKeys.size() * 16; i++)
		ttt += suc20BitKeys[i];
	vector<uint32_t> optional20BitKeys;
	optional20BitKeys.resize(ttt);
	ind = 0;
	for (uint32_t i = 0; i < optional16BitKeys.size(); i++)
		for (uint32_t j = 0; j < 16; j++)
			if (suc20BitKeys[16 * i + j])
				optional20BitKeys[ind++] = optional16BitKeys[i] | (j << 16);

	vector<bool> suc24BitKeys;
	suc24BitKeys.resize(ttt * 16);
	for (int i = 0; i < ttt * 16; i++)
		suc24BitKeys[i] = true;

	bool sbDone4 = false;
	bool sbDiff4 = false;
	bool sbDone5 = false;
	bool sbDiff5 = false;

	for (int kInd = 0; kInd < ttt * 16; kInd++)
	{
		//Create the appropriate key
		uint32_t ke = optional20BitKeys[kInd / 16];
		ke |= (((uint32_t)(kInd % 16)) << 20);

		for (int cind = 0; cind < NGOODPAIRS; cind++)
		{
			if (((c1s[cind] & hardMask) == (c2s[cind] & hardMask)) &&
				((c1s[cind] & ((uint64_t)(0x7800000000))) ^ (c2s[cind] & ((uint64_t)(0x7800000000)))))
			{
				uint32_t val1 = (ke&((uint32_t)0xffffff)) + (c1s[cind] & ((uint32_t)0xffffff));
				uint32_t val2 = (ke&((uint32_t)0xffffff)) + (c2s[cind] & ((uint32_t)0xffffff));
				if ((val1 >= 0x1000000) && (val2 >= 0x1000000))
				{
					suc24BitKeys[kInd] = 0;
					break;
				}
				if ((val1 < 0x1000000) && (val2 < 0x1000000))
				{
					suc24BitKeys[kInd] = 0;
					break;
				}
			}
		}
	}

	int NRemKeys = 0;
	for (int ll = 0; ll < ttt * 16; ll++)
		NRemKeys += suc24BitKeys[ll];
	printf("The number of remain keys after the first elimination is: %d\n", NRemKeys);

	for (int kInd = 0; kInd < ttt * 16; kInd++)
	{
		if (!suc24BitKeys[kInd])
			continue;
		uint32_t ke = optional20BitKeys[kInd / 16];
		ke |= (((uint32_t)(kInd % 16)) << 20);
		int inValsOutDiffs4[NGOODPAIRS][4] = {};
		int inValsOutDiffs5[NGOODPAIRS][4] = {};

		//Check consistency
		buildInValOutDiffs(4, c1s, c2s, ke, inValsOutDiffs4, 0, 0, 0);

		int inds4[17] = {};
		if (isWrongKey(inValsOutDiffs4, inds4, 0))
		{
			suc24BitKeys[kInd] = 0;
			continue;
		}

		buildInValOutDiffs(5, c1s, c2s, ke, inValsOutDiffs5, 0, 0, 0);
		int inds5[17] = {};
		if (isWrongKey(inValsOutDiffs5, inds5, 0))
		{
			suc24BitKeys[kInd] = 0;
			continue;
		}

		int rSboxTemp[8][16] = {};
		buildSbox(rSboxTemp, inValsOutDiffs4, inds4, 4, 0);
		buildSbox(rSboxTemp, inValsOutDiffs5, inds5, 5, 0);
		uint64_t stemp = 0;
		for (int val = 0; val < 16; val++)
			stemp |= ((uint64_t)rSboxTemp[4][val] << (4 * val));
		s4.push_back(stemp);
		stemp = 0;
		for (int val = 0; val < 16; val++)
			stemp |= ((uint64_t)rSboxTemp[5][val] << (4 * val));
		s5.push_back(stemp);
		resultKeys.push_back(ke);
		sbDone4 = 1;
		sbDone5 = 1;
		if (sbDone4 && !sbDiff4)
		{
			for (int i = 0; i < 16; i++)
			{
				if (((s4[0] & (((uint64_t)0xf) << 4 * i)) >> 4 * i) ^ rSboxTemp[4][i])
				{
					sbDiff4 = true;
					DS4++;
					break;
				}
			}
		}

		if (sbDone5 && !sbDiff5)
		{
			for (int i = 0; i < 16; i++)
			{
				if (((s5[0] & (((uint64_t)0xf) << 4 * i)) >> 4 * i) ^ rSboxTemp[5][i])
				{
					sbDiff5 = true;
					DS5++;
					break;
				}
			}
		}

		if (ke == theKey)
		{
			for (int exp = 4; exp < 6; exp++)
				for (int j = 1; j < 16; j++)
					if (rSboxTemp[exp][j] ^ rSboxTemp[exp][0] ^ s_block[exp][j] ^ s_block[exp][0])
					{
						printf("The first stage failed.\n");
						return false;
					}
			SUC1++;
		}
	}
	printf("The number of remain keys (out of 2^{24}) after the first stage is: %d\n", resultKeys.size());
	return true;
}

double choosediffBits(int diffBits[2], int s[16])
{
	int ddt[16][16] = {};
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt[i ^ j][s[i] ^ s[j]]++;
	double sump[4][4] = {};

	sump[0][0] = ((((double)ddt[1][1] / 16.)*0.5) + (((double)ddt[3][1] / 16.)*0.25) + (((double)ddt[7][1] / 16.)*0.125) + (((double)ddt[15][1] / 16.)*0.0625))
		*((((double)ddt[1][1] / 16.)*0.5) + (((double)ddt[3][1] / 16.)*0.25) + (((double)ddt[7][1] / 16.)*0.125) + (((double)ddt[15][1] / 16.)*0.0625))*0.5;
	sump[0][0] += ((((double)ddt[1][3] / 16.)*0.5) + (((double)ddt[3][3] / 16.)*0.25) + (((double)ddt[7][3] / 16.)*0.125) + (((double)ddt[15][3] / 16.)*0.0625))
		*((((double)ddt[1][3] / 16.)*0.5) + (((double)ddt[3][3] / 16.)*0.25) + (((double)ddt[7][3] / 16.)*0.125) + (((double)ddt[15][3] / 16.)*0.0625))*0.25;
	sump[0][0] += ((((double)ddt[1][7] / 16.)*0.5) + (((double)ddt[3][7] / 16.)*0.25) + (((double)ddt[7][7] / 16.)*0.125) + (((double)ddt[15][7] / 16.)*0.0625))
		*((((double)ddt[1][7] / 16.)*0.5) + (((double)ddt[3][7] / 16.)*0.25) + (((double)ddt[7][7] / 16.)*0.125) + (((double)ddt[15][7] / 16.)*0.0625))*0.125;
	sump[0][0] += ((((double)ddt[1][15] / 16.)*0.5) + (((double)ddt[3][15] / 16.)*0.25) + (((double)ddt[7][15] / 16.)*0.125) + (((double)ddt[15][15] / 16.)*0.0625))
		*((((double)ddt[1][15] / 16.)*0.5) + (((double)ddt[3][15] / 16.)*0.25) + (((double)ddt[7][15] / 16.)*0.125) + (((double)ddt[15][15] / 16.)*0.0625))*0.0625;

	sump[0][1] = ((((double)ddt[1][2] / 16.)*0.5) + (((double)ddt[3][2] / 16.)*0.25) + (((double)ddt[7][2] / 16.)*0.125) + (((double)ddt[15][2] / 16.)*0.0625))
		*((((double)ddt[1][2] / 16.)*0.5) + (((double)ddt[3][2] / 16.)*0.25) + (((double)ddt[7][2] / 16.)*0.125) + (((double)ddt[15][2] / 16.)*0.0625))*0.5;
	sump[0][1] += ((((double)ddt[1][6] / 16.)*0.5) + (((double)ddt[3][6] / 16.)*0.25) + (((double)ddt[7][6] / 16.)*0.125) + (((double)ddt[15][6] / 16.)*0.0625))
		*((((double)ddt[1][6] / 16.)*0.5) + (((double)ddt[3][6] / 16.)*0.25) + (((double)ddt[7][6] / 16.)*0.125) + (((double)ddt[15][6] / 16.)*0.0625))*0.25;
	sump[0][1] += ((((double)ddt[1][14] / 16.)*0.5) + (((double)ddt[3][14] / 16.)*0.25) + (((double)ddt[7][14] / 16.)*0.125) + (((double)ddt[15][14] / 16.)*0.0625))
		*((((double)ddt[1][14] / 16.)*0.5) + (((double)ddt[3][14] / 16.)*0.25) + (((double)ddt[7][14] / 16.)*0.125) + (((double)ddt[15][14] / 16.)*0.0625))*0.125;

	sump[0][2] = ((((double)ddt[1][4] / 16.)*0.5) + (((double)ddt[3][4] / 16.)*0.25) + (((double)ddt[7][4] / 16.)*0.125) + (((double)ddt[15][4] / 16.)*0.0625))
		*((((double)ddt[1][4] / 16.)*0.5) + (((double)ddt[3][4] / 16.)*0.25) + (((double)ddt[7][4] / 16.)*0.125) + (((double)ddt[15][4] / 16.)*0.0625))*0.5;
	sump[0][2] += ((((double)ddt[1][12] / 16.)*0.5) + (((double)ddt[3][12] / 16.)*0.25) + (((double)ddt[7][12] / 16.)*0.125) + (((double)ddt[15][12] / 16.)*0.0625))
		*((((double)ddt[1][12] / 16.)*0.5) + (((double)ddt[3][12] / 16.)*0.25) + (((double)ddt[7][12] / 16.)*0.125) + (((double)ddt[15][12] / 16.)*0.0625))*0.25;

	sump[0][3] = ((((double)ddt[1][8] / 16.)*0.5) + (((double)ddt[3][8] / 16.)*0.25) + (((double)ddt[7][8] / 16.)*0.125) + (((double)ddt[15][8] / 16.)*0.0625))
		*((((double)ddt[1][8] / 16.)*0.5) + (((double)ddt[3][8] / 16.)*0.25) + (((double)ddt[7][8] / 16.)*0.125) + (((double)ddt[15][8] / 16.)*0.0625))*0.5;

	sump[1][0] = ((((double)ddt[2][1] / 16.)*0.5) + (((double)ddt[6][1] / 16.)*0.25) + (((double)ddt[14][1] / 16.)*0.125))
		*((((double)ddt[2][1] / 16.)*0.5) + (((double)ddt[6][1] / 16.)*0.25) + (((double)ddt[14][1] / 16.)*0.125))*0.5;
	sump[1][0] += ((((double)ddt[2][3] / 16.)*0.5) + (((double)ddt[6][3] / 16.)*0.25) + (((double)ddt[14][3] / 16.)*0.125))
		*((((double)ddt[2][3] / 16.)*0.5) + (((double)ddt[6][3] / 16.)*0.25) + (((double)ddt[14][3] / 16.)*0.125))*0.25;
	sump[1][0] += ((((double)ddt[2][7] / 16.)*0.5) + (((double)ddt[6][7] / 16.)*0.25) + (((double)ddt[14][7] / 16.)*0.125))
		*((((double)ddt[2][7] / 16.)*0.5) + (((double)ddt[6][7] / 16.)*0.25) + (((double)ddt[14][7] / 16.)*0.125))*0.125;
	sump[1][0] += ((((double)ddt[2][15] / 16.)*0.5) + (((double)ddt[6][15] / 16.)*0.25) + (((double)ddt[14][15] / 16.)*0.125))
		*((((double)ddt[2][15] / 16.)*0.5) + (((double)ddt[6][15] / 16.)*0.25) + (((double)ddt[14][15] / 16.)*0.125))*0.0625;

	sump[1][1] = ((((double)ddt[2][2] / 16.)*0.5) + (((double)ddt[6][2] / 16.)*0.25) + (((double)ddt[14][2] / 16.)*0.125))
		*((((double)ddt[2][2] / 16.)*0.5) + (((double)ddt[6][2] / 16.)*0.25) + (((double)ddt[14][2] / 16.)*0.125))*0.5;
	sump[1][1] += ((((double)ddt[2][6] / 16.)*0.5) + (((double)ddt[6][6] / 16.)*0.25) + (((double)ddt[14][6] / 16.)*0.125))
		*((((double)ddt[2][6] / 16.)*0.5) + (((double)ddt[6][6] / 16.)*0.25) + (((double)ddt[14][6] / 16.)*0.125))*0.25;
	sump[1][1] += ((((double)ddt[2][14] / 16.)*0.5) + (((double)ddt[6][14] / 16.)*0.25) + (((double)ddt[14][14] / 16.)*0.125))
		*((((double)ddt[2][14] / 16.)*0.5) + (((double)ddt[6][14] / 16.)*0.25) + (((double)ddt[14][14] / 16.)*0.125))*0.125;

	sump[1][2] = ((((double)ddt[2][4] / 16.)*0.5) + (((double)ddt[6][4] / 16.)*0.25) + (((double)ddt[14][4] / 16.)*0.125))
		*((((double)ddt[2][4] / 16.)*0.5) + (((double)ddt[6][4] / 16.)*0.25) + (((double)ddt[14][4] / 16.)*0.125))*0.5;
	sump[1][2] += ((((double)ddt[2][12] / 16.)*0.5) + (((double)ddt[6][12] / 16.)*0.25) + (((double)ddt[14][12] / 16.)*0.125))
		*((((double)ddt[2][12] / 16.)*0.5) + (((double)ddt[6][12] / 16.)*0.25) + (((double)ddt[14][12] / 16.)*0.125))*0.25;

	sump[1][3] = ((((double)ddt[2][8] / 16.)*0.5) + (((double)ddt[6][8] / 16.)*0.25) + (((double)ddt[14][8] / 16.)*0.125))
		*((((double)ddt[2][8] / 16.)*0.5) + (((double)ddt[6][8] / 16.)*0.25) + (((double)ddt[14][8] / 16.)*0.125))*0.5;

	sump[2][0] = ((((double)ddt[4][1] / 16.)*0.5) + (((double)ddt[12][1] / 16.)*0.25))
		*((((double)ddt[4][1] / 16.)*0.5) + (((double)ddt[12][1] / 16.)*0.25))*0.5;
	sump[2][0] += ((((double)ddt[4][3] / 16.)*0.5) + (((double)ddt[12][3] / 16.)*0.25))
		*((((double)ddt[4][3] / 16.)*0.5) + (((double)ddt[12][3] / 16.)*0.25))*0.25;
	sump[2][0] += ((((double)ddt[4][7] / 16.)*0.5) + (((double)ddt[12][7] / 16.)*0.25))
		*((((double)ddt[4][7] / 16.)*0.5) + (((double)ddt[12][7] / 16.)*0.25))*0.125;
	sump[2][0] += ((((double)ddt[4][15] / 16.)*0.5) + (((double)ddt[12][15] / 16.)*0.25))
		*((((double)ddt[4][15] / 16.)*0.5) + (((double)ddt[12][15] / 16.)*0.25))*0.0625;

	sump[2][1] = ((((double)ddt[4][2] / 16.)*0.5) + (((double)ddt[12][2] / 16.)*0.25))
		*((((double)ddt[4][2] / 16.)*0.5) + (((double)ddt[12][2] / 16.)*0.25))*0.5;
	sump[2][1] += ((((double)ddt[4][6] / 16.)*0.5) + (((double)ddt[12][6] / 16.)*0.25))
		*((((double)ddt[4][6] / 16.)*0.5) + (((double)ddt[12][6] / 16.)*0.25))*0.25;
	sump[2][1] += ((((double)ddt[4][14] / 16.)*0.5) + (((double)ddt[12][14] / 16.)*0.25))
		*((((double)ddt[4][14] / 16.)*0.5) + (((double)ddt[12][14] / 16.)*0.25))*0.125;

	sump[2][2] = ((((double)ddt[4][4] / 16.)*0.5) + (((double)ddt[12][4] / 16.)*0.25))
		*((((double)ddt[4][4] / 16.)*0.5) + (((double)ddt[12][4] / 16.)*0.25))*0.5;
	sump[2][2] += ((((double)ddt[4][12] / 16.)*0.5) + (((double)ddt[12][12] / 16.)*0.25))
		*((((double)ddt[4][12] / 16.)*0.5) + (((double)ddt[12][12] / 16.)*0.25))*0.25;

	sump[2][3] = ((((double)ddt[4][8] / 16.)*0.5) + (((double)ddt[12][8] / 16.)*0.25))
		*((((double)ddt[4][8] / 16.)*0.5) + (((double)ddt[12][8] / 16.)*0.25))*0.5;


	sump[3][0] = ((((double)ddt[8][1] / 16.)*0.5))
		*((((double)ddt[8][1] / 16.)*0.5))*0.5;
	sump[3][0] += ((((double)ddt[8][3] / 16.)*0.5))
		*((((double)ddt[8][3] / 16.)*0.5))*0.25;
	sump[3][0] += ((((double)ddt[8][7] / 16.)*0.5))
		*((((double)ddt[8][7] / 16.)*0.5))*0.125;
	sump[3][0] += ((((double)ddt[8][15] / 16.)*0.5))
		*((((double)ddt[8][15] / 16.)*0.5))*0.0625;

	sump[3][1] = ((((double)ddt[8][2] / 16.)*0.5))
		*((((double)ddt[8][2] / 16.)*0.5))*0.5;
	sump[3][1] += ((((double)ddt[8][6] / 16.)*0.5))
		*((((double)ddt[8][6] / 16.)*0.5))*0.25;
	sump[3][1] += ((((double)ddt[8][14] / 16.)*0.5))
		*((((double)ddt[8][14] / 16.)*0.5))*0.125;

	sump[3][2] = ((((double)ddt[8][4] / 16.)*0.5))
		*((((double)ddt[8][4] / 16.)*0.5))*0.5;
	sump[3][2] += ((((double)ddt[8][12] / 16.)*0.5))
		*((((double)ddt[8][12] / 16.)*0.5))*0.25;

	sump[3][3] = ((((double)ddt[8][8] / 16.)*0.5))
		*((((double)ddt[8][8] / 16.)*0.5))*0.5;

	double maxp = sump[0][0];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (sump[i][j] > maxp)
			{

				maxp = sump[i][j];
				diffBits[0] = i;
				diffBits[1] = j;
			}
		}
	}

	return maxp;
}

bool fndAddS(int db2[2], vector<uint32_t> &resultKeys, uint32_t the24BitsKey,
	uint32_t key1[8], vector<uint64_t> &s1, vector<uint64_t> &s2, vector<uint64_t> &s4,
	vector<uint64_t> &s5, uint64_t c1s[NGOODPAIRS], uint64_t c2s[NGOODPAIRS])
{
	NEcr = 0;
	int NB[25] = { 44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,63,4,5,6,7,8,9,10,11,12 };

	if (!s4.size())
		return 0;
	int s_4[16] = {};
	uint64_t stemp = s4[0];
	for (int val = 0; val < 16; val++)
	{
		s_4[val] = stemp % 16;
		stemp >>= 4;
	}
	choosediffBits(db2, s_4);

	uint32_t key2[8] = {};
	for (int i = 0; i < 8; i++)
		key2[i] = key1[i];
	key2[0] ^= ((uint32_t)1 << 16 + db2[0]);
	key2[1] ^= ((uint32_t)1 << 27 + db2[1]);
	key2[2] ^= ((uint32_t)1 << 16 + db2[0]);

	uint64_t hardMask = 0x000007f8fffff007;
	uint64_t weakMask = 0x000007f8ffff8007;
	bool notFound = true;

	//cout << "resK = " << resultKeys.size() << endl;
	while ((notFound) && (NEcr < 0x8000000))
	{
		uint64_t p1 = 0;
		for (int i = 0; i < 4; i++)
			p1 |= ((uint64_t)(rand()) << (15 * i));
		p1 |= ((uint64_t)(rand() % 16) << 60);

		uint64_t c1 = encrypt(p1, key1, s_block, 32);
		uint64_t c2 = encrypt(p1, key2, s_block, 32);

		uint32_t numGoodNewPln = 0;

		if ((c1&weakMask) == (c2&weakMask))
		{
			c1s[numGoodNewPln] = c1;
			c2s[numGoodNewPln++] = c2;
			notFound = false;

			uint32_t plnInd = 1;
			while (numGoodNewPln < NGOODPAIRS)
			{
				//if (NEcr > 0x8000000)
				if ((numGoodNewPln == 1 && plnInd == 0x40000) || (plnInd > 0x8000000))
				{
					notFound = true;
					break;
				}
				uint64_t newP1 = p1;
				uint32_t temp = plnInd++;
				int NBInd = 0;
				while (temp)
				{
					newP1 ^= ((uint64_t)(temp % 2) << NB[NBInd++]);
					temp >>= 1;
				}
				c1 = encrypt(newP1, key1, s_block, 32);
				c2 = encrypt(newP1, key2, s_block, 32);

				if ((c1&weakMask) == (c2&weakMask))
				{
					c1s[numGoodNewPln] = c1;
					c2s[numGoodNewPln++] = c2;
					int nwk1 = 0;
					for (uint32_t k = 0; k < resultKeys.size(); k++)
					{
						//uint32_t ke = resultKeys[k];
						//Difference in S_1 -> there is carry to bit 23 in single addition.
						if (((c1&hardMask) == (c2&hardMask)) &&
							((c1&((uint64_t)(0x780000000000000))) ^ (c2&((uint64_t)(0x780000000000000)))))
						{
							uint32_t val1 = (resultKeys[k] & (uint32_t)0xfff) + (c1&((uint32_t)0xfff));
							uint32_t val2 = (resultKeys[k] & (uint32_t)0xfff) + (c2&((uint32_t)0xfff));
							if ((val1 >= 4096) && (val2 >= 4096))
							{
								resultKeys.erase(resultKeys.begin() + k);
								s4.erase(s4.begin() + k);
								s5.erase(s5.begin() + (k--));
								nwk1++;
								continue;
							}
							if ((val1 < 4096) && (val2 < 4096))
							{
								resultKeys.erase(resultKeys.begin() + k);
								s4.erase(s4.begin() + k);
								s5.erase(s5.begin() + (k--));								nwk1++;
								continue;
							}
						}
					}
				}
			}
			if (notFound)
				break;
		}
	}

	if (notFound)
	{
		printf("Not found good pairs for the second stage.\n");
		return false;
	}
	printf("The number of remain keys after the advance elimination of the second stage is: %d\n", resultKeys.size());

	bool sbDone2 = false;
	bool sbDiff2 = false;
	bool sbDone1 = false;
	bool sbDiff1 = false;

	for (int i = 0; i < resultKeys.size(); i++)
	{
		uint32_t ke = resultKeys[i];
		int inValsOutDiffs1[NGOODPAIRS][4] = {};
		int inValsOutDiffs2[NGOODPAIRS][4] = {};

		//Check consistency
		buildInValOutDiffs(1, c1s, c2s, ke, inValsOutDiffs1, 0, 0, 0);
		int inds1[17] = {};
		if (isWrongKey(inValsOutDiffs1, inds1, 0))
		{
			resultKeys.erase(resultKeys.begin() + i);
			s4.erase(s4.begin() + i);
			s5.erase(s5.begin() + (i--));
			continue;
		}

		buildInValOutDiffs(2, c1s, c2s, ke, inValsOutDiffs2, 0, 0, 0);
		int inds2[17] = {};
		if (isWrongKey(inValsOutDiffs2, inds2, 0))
		{
			resultKeys.erase(resultKeys.begin() + i);
			s4.erase(s4.begin() + i);
			s5.erase(s5.begin() + (i--));
			continue;
		}

		int rSboxTemp[8][16] = {};
		buildSbox(rSboxTemp, inValsOutDiffs1, inds1, 1, 0);
		buildSbox(rSboxTemp, inValsOutDiffs2, inds2, 2, 0);
		stemp = 0;
		for (int val = 0; val < 16; val++)
			stemp |= ((uint64_t)rSboxTemp[1][val] << (4 * val));
		s1.push_back(stemp);
		stemp = 0;
		for (int val = 0; val < 16; val++)
			stemp |= ((uint64_t)rSboxTemp[2][val] << (4 * val));
		s2.push_back(stemp);

		sbDone2 = 1;
		sbDone1 = 1;
		if (sbDone2 && !sbDiff2)
		{
			for (int j = 0; j < 16; j++)
			{
				if (((s2[0] & (((uint64_t)0xf) << 4 * j)) >> 4 * j) ^ rSboxTemp[2][j])
				{
					sbDiff2 = true;
					DS2++;
					break;
				}
			}
		}

		if (sbDone1 && !sbDiff1)
		{
			for (int j = 0; j < 16; j++)
			{
				if (((s1[0] & (((uint64_t)0xf) << 4 * j)) >> 4 * j) ^ rSboxTemp[1][j])
				{
					sbDiff1 = true;
					DS1++;
					break;
				}
			}
		}
		if (ke == the24BitsKey)
		{
			for (int sbox = 1; sbox < 3; sbox++)
				for (int j = 0; j < 16; j++)
					if (rSboxTemp[sbox][j] ^ s_block[sbox][j] ^ rSboxTemp[sbox][0] ^ s_block[sbox][0])
					{
						printf("The second stage failed.\n");
						return false;
					}
			SUC2++;
		}
	}
	printf("The number of remain keys after the second stage is: %d\n", resultKeys.size());
	return true;
}

bool createThirdPCSet(int db3[2], uint32_t key1[8], vector<uint64_t> s2, uint64_t c1s[NGOODPAIRS], uint64_t c2s[NGOODPAIRS], uint64_t ps[NGOODPAIRS])
{
	NEcr = 0;
	if (!s2.size())
		return 0;
	int NB[25] = { 28,29,30,31,0,1,2,3,4,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55 };

	bool goodDiffs[16] = {};

	int ddt[16][16] = {};
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
			ddt[i ^ j][s_block[6][i] ^ s_block[6][j]]++;

	int s_2[16] = {};
	uint64_t stemp = s2[0];
	for (int val = 0; val < 16; val++)
	{
		s_2[val] = stemp % 16;
		stemp >>= 4;
	}
	choosediffBits(db3, s_2);

	uint32_t key2[8] = {};
	for (int i = 0; i < 8; i++)
		key2[i] = key1[i];
	key2[0] ^= ((uint32_t)1 << 8 + db3[0]);
	key2[1] ^= ((uint32_t)1 << 19 + db3[1]);
	key2[2] ^= ((uint32_t)1 << 8 + db3[0]);

	uint64_t weakMask = 0xf800000707ffff80;
	bool notFound = true;
	//uint32_t numOutIter = 0;
	while ((notFound) && (NEcr < 0x8000000))
	{
		//numOutIter++;
		uint64_t p1 = 0;
		for (int i = 0; i < 4; i++)
			p1 |= ((uint64_t)(rand()) << (15 * i));
		p1 |= ((uint64_t)(rand() % 16) << 60);

		uint64_t c1 = encrypt(p1, key1, s_block, 32);
		uint64_t c2 = encrypt(p1, key2, s_block, 32);

		uint32_t numGoodNewPln = 0;

		if ((c1&weakMask) == (c2&weakMask))
		{
			ps[numGoodNewPln] = p1;
			c1s[numGoodNewPln] = c1;
			c2s[numGoodNewPln++] = c2;
			notFound = false;

			uint32_t plnInd = 1;
			while (numGoodNewPln < NGOODPAIRS)
			{
				//if (plnInd == 1)
				//	cout << "yoffi" << endl;

				//if (NEcr > 0x8000000)
				if ((numGoodNewPln == 1 && plnInd == 0x40000) || (plnInd > 0x8000000))
				{
					notFound = true;
					break;
				}
				uint64_t newP1 = p1;
				uint32_t temp = plnInd++;
				int NBInd = 0;
				while (temp)
				{
					newP1 ^= ((uint64_t)(temp % 2) << NB[NBInd++]);
					temp >>= 1;
				}
				c1 = encrypt(newP1, key1, s_block, 32);
				c2 = encrypt(newP1, key2, s_block, 32);

				if ((c1&weakMask) == (c2&weakMask))
				{
					ps[numGoodNewPln] = newP1;
					c1s[numGoodNewPln] = c1;
					c2s[numGoodNewPln++] = c2;
					int nwk1 = 0;
				}
			}
			if (notFound)
				break;
		}
	}

	if (notFound)
	{
		//file << "oish" << endl;
		return false;
	}
	return true;
}

bool lastElimination07(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS],
	short sb[8][16], uint32_t ke, bool LSB8_A[256], bool k2_8[256], int db)
{
	bool LSB8_Atemp[256] = {};
	for (int k2Op = 0; k2Op < 256; k2Op++)
	{
		for (int i = 0; i < 256; i++)
			LSB8_Atemp[i] = 1;
		for (int lsbOp = 0; lsbOp < 256; lsbOp++)
		{
			for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
			{
				int cxor = (((c1s_1[cInd] & 0x78000) ^ (c2s_1[cInd] & 0x78000)) >> 15);
				if (!cxor)
					continue;
				uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
				uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << 31), sb);
				uint32_t c1_31R = c1_31 & 0xffffffff;
				uint32_t c2_31R = c2_31 & 0xffffffff;
				c1_31R ^= lsbOp;
				c2_31R ^= lsbOp;
				c1_31R += k2Op;
				if (db)
					c2_31R += k2Op;
				else
					c2_31R += (k2Op ^ 128);
				int xorCandidate = (sb[1][((c1_31R & 0xf0) >> 4)]) ^ (sb[1][((c2_31R & 0xf0) >> 4)]);
				if (cxor^xorCandidate)
				{
					LSB8_Atemp[lsbOp] = 0;
					break;
				}
			}
		}
		for (int i = 0; i < 256; i++)
			if (LSB8_Atemp[i])
			{
				k2_8[k2Op] = 1;
				LSB8_A[i] = 1;
			}
	}
	for (int k = 0; k < 256; k++)
		if (k2_8[k])
			return true;
	return false;
}

bool lastElimination811(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS], short sb[8][16],
	uint32_t ke, vector<bool> &LSB12_A, vector<bool> &k2_12, int db,
	vector <uint32_t> As, vector <uint32_t> k2s)
{
	vector<bool> LSB12_Atemp(LSB12_A.size());
	for (int k2Op = 0; k2Op < k2_12.size(); k2Op++)
	{
		k2_12[k2Op] = 0;
		uint32_t k2Practice = (((k2Op % 16) << 8) | (k2s[k2Op / 16]));
		for (int i = 0; i < LSB12_A.size(); i++)
			LSB12_Atemp[i] = 1;
		for (int lsbOp = 0; lsbOp < LSB12_A.size(); lsbOp++)
		{
			uint32_t lsbPractice = (((lsbOp % 16) << 8) | (As[lsbOp / 16]));
			for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
			{
				int cxor = (((c1s_1[cInd] & 0x780000) ^ (c2s_1[cInd] & 0x780000)) >> 19);
				if (!cxor)
					continue;
				uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
				uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << 31), sb);
				uint32_t c1_31R = c1_31 & 0xffffffff;
				uint32_t c2_31R = c2_31 & 0xffffffff;

				c1_31R ^= lsbPractice;
				c2_31R ^= lsbPractice;
				c1_31R += k2Practice;
				c2_31R += (k2Practice ^ (128 << db));
				int xorCandidate = (sb[2][((c1_31R & 0xf00) >> 8)]) ^
					(sb[2][((c2_31R & 0xf00) >> 8)]);
				if (cxor^xorCandidate)
				{
					LSB12_Atemp[lsbOp] = 0;
					break;
				}
			}
		}
		for (int i = 0; i < LSB12_A.size(); i++)
			if (LSB12_Atemp[i])
			{
				k2_12[k2Op] = 1;
				LSB12_A[i] = 1;
			}
	}
	for (int k = 0; k < k2_12.size(); k++)
		if (k2_12[k])
			return true;
	return false;
}

bool lastElimination1219(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS], short sb[8][16],
	uint32_t ke, vector<bool> &LSB20_A, vector<bool> &k2_20, int db[2],
	vector <uint32_t> As, vector <uint32_t> k2s)
{
	vector<bool> LSB20_Atemp(LSB20_A.size());
	for (int k2Op = 0; k2Op < k2_20.size(); k2Op++)
	{
		k2_20[k2Op] = 0;
		uint32_t k2Practice = (((k2Op % 256) << 12) | (k2s[k2Op / 256]));
		for (int i = 0; i < LSB20_A.size(); i++)
			LSB20_Atemp[i] = 1;
		for (int lsbOp = 0; lsbOp < LSB20_A.size(); lsbOp++)
		{
			uint32_t lsbPractice = (((lsbOp % 256) << 12) | (As[lsbOp / 256]));
			for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
			{
				int cxor = (((c1s_1[cInd] & 0x78000000) ^ (c2s_1[cInd] & 0x78000000)) >> 27);
				if (!cxor)
					continue;
				uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
				uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << (8 + db[0])), sb);
				uint32_t c1_31R = c1_31 & 0xffffffff;
				uint32_t c2_31R = c2_31 & 0xffffffff;

				c1_31R ^= lsbPractice;
				c2_31R ^= lsbPractice;
				c1_31R += k2Practice;
				c2_31R += (k2Practice ^ ((uint32_t)1 << (19 + db[1])));
				int xorCandidate = (sb[4][((c1_31R & 0xf0000) >> 16)]) ^
					(sb[4][((c2_31R & 0xf0000) >> 16)]);
				if (cxor^xorCandidate)
				{
					LSB20_Atemp[lsbOp] = 0;
					break;
				}
			}
		}
		for (int i = 0; i < LSB20_A.size(); i++)
			if (LSB20_Atemp[i])
			{
				k2_20[k2Op] = 1;
				LSB20_A[i] = 1;
			}
	}
	for (int k = 0; k < k2_20.size(); k++)
		if (k2_20[k])
			return true;
	return false;
}

bool lastElimination2023(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS], short sb[8][16],
	uint32_t ke, vector<bool> &LSB24_A, vector<bool> &k2_24, int db[2],
	vector <uint32_t> As, vector <uint32_t> k2s)
{
	vector<bool> LSB24_Atemp(LSB24_A.size());
	for (int k2Op = 0; k2Op < k2_24.size(); k2Op++)
	{
		k2_24[k2Op] = 0;
		uint32_t k2Practice = (((k2Op % 16) << 20) | (k2s[k2Op / 16]));
		for (int i = 0; i < LSB24_A.size(); i++)
			LSB24_Atemp[i] = 1;
		for (int lsbOp = 0; lsbOp < LSB24_A.size(); lsbOp++)
		{
			uint32_t lsbPractice = (((lsbOp % 16) << 20) | (As[lsbOp / 16]));
			for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
			{
				int cxor = ((c1s_1[cInd] & 7) ^ (c2s_1[cInd] & 7));
				cxor <<= 1;
				if ((c1s_1[cInd] & 0x80000000) ^ (c2s_1[cInd] & 0x80000000))
					cxor |= 1;
				if (!cxor)
					continue;
				uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
				uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << (8 + db[0])), sb);
				uint32_t c1_31R = c1_31 & 0xffffffff;
				uint32_t c2_31R = c2_31 & 0xffffffff;

				c1_31R ^= lsbPractice;
				c2_31R ^= lsbPractice;
				c1_31R += k2Practice;
				c2_31R += (k2Practice ^ ((uint32_t)1 << (19 + db[1])));

				int xorCandidate = (sb[5][((c1_31R & 0xf00000) >> 20)]) ^
					(sb[5][((c2_31R & 0xf00000) >> 20)]);
				if (cxor^xorCandidate)
				{
					LSB24_Atemp[lsbOp] = 0;
					break;
				}
			}
		}
		for (int i = 0; i < LSB24_A.size(); i++)
			if (LSB24_Atemp[i])
			{
				k2_24[k2Op] = 1;
				LSB24_A[i] = 1;
			}
	}
	for (int k = 0; k < k2_24.size(); k++)
		if (k2_24[k])
			return true;
	return false;
}

bool lastElimination2431(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS], short sb[8][16],
	uint32_t ke, vector<bool> &LSB32_A, vector<bool> &k2_32, int db[2],
	vector <uint32_t> As, vector <uint32_t> k2s)
{
	vector<bool> LSB32_Atemp(LSB32_A.size());
	for (int k2Op = 0; k2Op < k2_32.size(); k2Op++)
	{
		k2_32[k2Op] = 0;
		uint32_t k2Practice = (((k2Op % 256) << 24) | (k2s[k2Op / 256]));
		for (int i = 0; i < LSB32_A.size(); i++)
			LSB32_Atemp[i] = 1;
		for (int lsbOp = 0; lsbOp < LSB32_A.size(); lsbOp++)
		{
			uint32_t lsbPractice = (((lsbOp % 256) << 24) | (As[lsbOp / 256]));
			for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
			{
				int cxor = ((c1s_1[cInd] & 0x7f8) ^ (c2s_1[cInd] & 0x7f8));
				if (cxor)
				{
					cxor >>= 3;
					uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
					uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << (16 + db[0])), sb);
					uint32_t c1_31R = c1_31 & 0xffffffff;
					uint32_t c2_31R = c2_31 & 0xffffffff;

					c1_31R ^= lsbPractice;
					c2_31R ^= lsbPractice;
					c1_31R += k2Practice;
					c2_31R += (k2Practice ^ ((uint32_t)1 << (27 + db[1])));

					int xorCandidate = (((sb[7][((c1_31R & 0xf0000000) >> 28)]) ^
						(sb[7][((c2_31R & 0xf0000000) >> 28)])) << 4);
					xorCandidate |= ((sb[6][((c1_31R & 0xf000000) >> 24)]) ^
						(sb[6][((c2_31R & 0xf000000) >> 24)]));
					if (cxor^xorCandidate)
					{
						LSB32_Atemp[lsbOp] = 0;
						break;
					}
				}
			}
		}
		for (int i = 0; i < LSB32_A.size(); i++)
			if (LSB32_Atemp[i])
			{
				k2_32[k2Op] = 1;
				LSB32_A[i] = 1;
			}
	}
	for (int k = 0; k < k2_32.size(); k++)
		if (k2_32[k])
			return true;
	return false;
}

int lastEliminationR30(uint64_t c1s_1[NGOODPAIRS], uint64_t c2s_1[NGOODPAIRS],
	uint64_t c1s_2[NGOODPAIRS], uint64_t c2s_2[NGOODPAIRS], uint64_t c1s_3[NGOODPAIRS], uint64_t c2s_3[NGOODPAIRS],
	short sb[8][16], uint32_t ke, int db1, int db2[2], int db3[2],
	vector <uint32_t> &As, vector <uint32_t> &k2s, vector<uint64_t> &full_k2_A, vector<uint32_t> &k3Canddt,
	uint32_t k1, uint32_t k3)
{
	bool k3Temp[4096] = {};
	for (int k2Ind = 0; k2Ind < k2s.size(); k2Ind++)
	{
		for (int AInd = 0; AInd < As.size(); AInd++)
		{
			for (int i = 0; i < 4096; i++)
				k3Temp[i] = 1;
			for (int k3Ind = 0; k3Ind < 4096; k3Ind++)
			{
				for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
				{
					uint64_t c1_31 = round(c1s_3[cInd], ke, sb);
					c1_31 ^= As[AInd];
					uint64_t c2_31 = round(c2s_3[cInd], ke ^ ((uint32_t)1 << (8 + db3[0])), sb);
					c2_31 ^= As[AInd];
					uint64_t c1_30 = round(c1_31, k2s[k2Ind], sb);
					c1_30 ^= As[AInd];
					uint64_t c2_30 = round(c2_31, k2s[k2Ind] ^ ((uint32_t)1 << (19 + db3[1])), sb);
					c2_30 ^= As[AInd];
					//cout << c1_30 << "\t" << c2_30 << endl;
					int cxor = (((c1_30 & (uint64_t)0x78000000000000) ^ (c2_30 & (uint64_t)0x78000000000000)) >> 51);
					//cout << cxor << endl;
					if (cxor)
					{
						c1_30 += k3Ind;
						c2_30 += (k3Ind ^ ((uint32_t)1 << (8 + db3[0])));

						int xorCandidate = ((sb[2][((c1_30 & 0xf00) >> 8)]) ^
							(sb[2][((c2_30 & 0xf00) >> 8)]));
						//xorCandidate |= ((((sb[4][((c1_30 & 0xf000) >> 12)]) ^
						//	(sb[4][((c2_30 & 0xf000) >> 12)]))) << 4);
						//if ((ke == k1) && (k3Ind == (k3 & 0xfff)))
						//	cout << cxor << "\t" << xorCandidate << endl;
						if (cxor^xorCandidate)
						{
							k3Temp[k3Ind] = 0;
							break;
						}
					}
				}
			}
			vector<uint32_t> k3tempVctr0;
			for (int i = 0; i < 4096; i++)
				if (k3Temp[i])
					k3tempVctr0.push_back(i);
			//cout << k3Canddt.size() << endl;
			vector<uint32_t> k3tempVctr1;
			for (uint32_t k3Ind = 0; k3Ind < 256; k3Ind++)
			{
				for (int k3Old = 0; k3Old < k3tempVctr0.size(); k3Old++)
				{
					bool goodCanddt = 1;
					for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
					{
						uint64_t c1_31 = round(c1s_2[cInd], ke, sb);
						c1_31 ^= As[AInd];
						uint64_t c2_31 = round(c2s_2[cInd], ke ^ ((uint32_t)1 << (16 + db2[0])), sb);
						c2_31 ^= As[AInd];
						uint64_t c1_30 = round(c1_31, k2s[k2Ind], sb);
						c1_30 ^= As[AInd];
						uint64_t c2_30 = round(c2_31, k2s[k2Ind] ^ ((uint32_t)1 << (27 + db2[1])), sb);
						c2_30 ^= As[AInd];
						//cout << c1_30 << "\t" << c2_30 << endl;
						uint32_t cxor = (((c1_30 & (uint64_t)0x7800000000000000) ^ (c2_30 & (uint64_t)0x7800000000000000)) >> 59);
						//cxor |= (((c1_30 & (uint64_t)0x700000000) ^ (c2_30 & (uint64_t)0x700000000)) >> 27);
						if (cxor)
						{
							//uint32_t c1_30R = c1_30 & 0xffffffff;
							//uint32_t c2_30R = c2_30 & 0xffffffff;
							c1_30 += ((k3Ind << 12) | (k3tempVctr0[k3Old]));
							c2_30 += (((k3Ind << 12) | (k3tempVctr0[k3Old])) ^ ((uint32_t)1 << (16 + db2[0])));

							int xorCandidate = ((sb[4][((c1_30 & 0xf0000) >> 16)]) ^
								(sb[4][((c2_30 & 0xf0000) >> 16)]));
							//xorCandidate |= ((((sb[2][((c1_30 & 0xf00000) >> 20)]) ^
							//	(sb[2][((c2_30 & 0xf00000) >> 20)]))) << 4);
							if (cxor^xorCandidate)
							{
								goodCanddt = 0;
								break;
							}
						}
					}
					if (goodCanddt)
						k3tempVctr1.push_back(((k3Ind << 12) | (k3tempVctr0[k3Old])));
				}
			}
			//cout << k3Temp.size() << endl;
			k3tempVctr0.clear();
			for (int k3Ind = 0; k3Ind < 4096; k3Ind++)
			{
				for (int k3Old = 0; k3Old < k3tempVctr1.size(); k3Old++)
				{
					bool goodCanddt = 1;
					for (int cInd = 0; cInd < NGOODPAIRS; cInd++)
					{
						uint64_t c1_31 = round(c1s_1[cInd], ke, sb);
						c1_31 ^= As[AInd];
						uint64_t c2_31 = round(c2s_1[cInd], ke ^ ((uint32_t)1 << 31), sb);
						c2_31 ^= As[AInd];
						uint64_t c1_30 = round(c1_31, k2s[k2Ind], sb);
						c1_30 ^= As[AInd];
						uint64_t c2_30 = round(c2_31, k2s[k2Ind] ^ ((uint32_t)1 << (7 + db1)), sb);
						c2_30 ^= As[AInd];
						//cout << c1_30 << "\t" << c2_30 << endl;
						int cxor = (((c1_30 & (uint64_t)0x78000000000) ^ (c2_30 & (uint64_t)0x78000000000)) >> 39);
						//cout << cxor << endl;
						if (cxor)
						{
							c1_30 += ((k3Ind << 20) | (k3tempVctr1[k3Old]));
							c2_30 += (((k3Ind << 20) | (k3tempVctr1[k3Old])) ^ ((uint32_t)1 << 31));

							int xorCandidate = ((sb[7][((c1_30 & 0xf0000000) >> 28)]) ^
								(sb[7][((c2_30 & 0xf0000000) >> 28)]));
							if (cxor^xorCandidate)
							{
								goodCanddt = 0;
								break;
							}
						}
					}
					if (goodCanddt)
					{
						k3tempVctr0.push_back(((k3Ind << 20) | (k3tempVctr1[k3Old])));
						bool isin = false;
						uint32_t val = ((k3Ind << 20) | (k3tempVctr1[k3Old]));
						for (int i = 0; i < k3Canddt.size(); i++)
							if (k3Canddt[i] == val)
								isin = 1;
						if (!isin)
							k3Canddt.push_back(((k3Ind << 20) | (k3tempVctr1[k3Old])));
					}
				}
			}
			//cout << k3Canddt.size() << endl;
			if (k3tempVctr0.size())
				full_k2_A.push_back((((uint64_t)k2s[k2Ind]) << 32) | ((uint64_t)As[AInd]));
		}
	}
	return full_k2_A.size();
}

void keySboxRecoveryAttackZD()
{
	int NB[30] = { 52,53,54,55,56,57,58,59,60,61,62,63,12,13,14,15,16,17,18,19,20,32,33,34,35,36,37,38,21,22 };
	uint64_t hardMask = 0x7ff800ff807fff;
	uint64_t harderMask = 0x7ff800fff07fff;
	uint64_t weakMask = 0x7ff80080007fff;
	DS5 = 0;
	DS4 = 0;
	DS1 = 0;
	DS2 = 0;

	for (int expInd = 0; expInd < 100; expInd++)
	{
		int NremKbfrLastElim = 0;
		NEcr = 0;
		printf("Experiment num. %d\n", expInd);
		printf("Time = %d\n", time(0));
		for (int s = 0; s < 8; s++)
		{
			int tempArr[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
			randomShuffle(tempArr, 16);
			for (int val = 0; val < 16; val++)
				s_block[s][val] = tempArr[val];
		}
		printf("The generated S-boxes are:\n");
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 16; j++)
				printf("%x\t", s_block[i][j]);
			printf("\n");
		}


		uint64_t c1s_1[NGOODPAIRS] = {};
		uint64_t c2s_1[NGOODPAIRS] = {};
		uint64_t c1s_2[NGOODPAIRS] = {};
		uint64_t c2s_2[NGOODPAIRS] = {};
		uint64_t c1s_3[NGOODPAIRS] = {};
		uint64_t c2s_3[NGOODPAIRS] = {};
		uint64_t ps[NGOODPAIRS] = {};

		uint32_t key1[8] = {};
		for (int i = 0; i < 8; i++)
			key1[i] = ((uint64_t)rand()) | ((uint64_t)(rand()) << 15) | ((uint64_t)(rand() % 4) << 30);
		printf("The generated sub-keys K1 and K2 are: %lx\t%lx\n", key1[0], key1[1]);

		uint32_t the19BitsKey = ((key1[0] & ((uint32_t)0x7ffff0)) >> 4);
		uint32_t the24BitsKey = (key1[0] & ((uint32_t)0xffffff));
		uint32_t theEntireKey = key1[0];

		bool suc = false;
		int DB = 0;

		bool notFound = true;
		//Find good pairs for the first stage.
		while ((notFound) && (NEcr < 0x8000000))
		{
			uint64_t p1 = 0;
			for (int i = 0; i < 4; i++)
				p1 |= ((uint64_t)(rand()) << (15 * i));
			p1 |= ((uint64_t)(rand() % 16) << 60);

			uint64_t c1 = encrypt(p1, key1, s_block, 32);
			uint64_t c2[4] = {};

			for (int diffBit = 0; diffBit < 4; diffBit++)
			{
				if (suc)
					break;

				uint32_t key2[8] = {};
				for (int i = 0; i < 8; i++)
					key2[i] = key1[i];
				key2[0] ^= ((uint32_t)1 << 31);
				key2[1] ^= ((uint32_t)1 << (7 + diffBit));
				key2[2] ^= ((uint32_t)1 << 31);

				c2[diffBit] = encrypt(p1, key2, s_block, 32);

				uint32_t numGoodNewPln = 0;

				if ((c1&weakMask) == (c2[diffBit] & weakMask))
				{
					DB = diffBit;
					c1s_1[numGoodNewPln] = c1;
					c2s_1[numGoodNewPln++] = c2[diffBit];
					notFound = false;
					suc = true;

					uint32_t plnInd = 1;
					while (numGoodNewPln < NGOODPAIRS)
					{
						if ((numGoodNewPln == 1 && plnInd == 0x800000) || (plnInd > 0x8000000))
						{
							notFound = true;
							suc = false;
							break;
						}
						uint64_t newP1 = p1;
						uint32_t temp = plnInd++;
						int NBInd = 0;
						while (temp)
						{
							newP1 ^= ((uint64_t)(temp % 2) << NB[NBInd++]);
							temp >>= 1;
						}
						c1 = encrypt(newP1, key1, s_block, 32);
						uint64_t newC2 = encrypt(newP1, key2, s_block, 32);

						if ((c1&weakMask) == (newC2&weakMask))
						{
							c1s_1[numGoodNewPln] = c1;
							c2s_1[numGoodNewPln++] = newC2;
						}
					}
					if (!suc)
						continue;
					if (notFound)
						break;
				}
			}
		}

		if (notFound)
		{
			printf("Failure: not enough pairs found.\n");
			continue;
		}

		vector<uint64_t> s1(0);
		vector<uint64_t> s2(0);
		vector<uint64_t> s4(0);
		vector<uint64_t> s5(0);
		vector<uint32_t> resultKeys(0);
		int ddt[16][16] = {};
		for (int i = 0; i < 16; i++)
			for (int j = 0; j < 16; j++)
				ddt[i ^ j][s_block[0][i] ^ s_block[0][j]]++;
		//The first stage.
		if (!SboxRecoveryAttack0(c1s_1, c2s_1,
			the24BitsKey, resultKeys, s4, s5))
		{
			printf("The first stage is failed.\n");
			continue;
		}

		printf("The number of GOST encriptions used in the first stage is: %d\n", NEcr);

		//The second stage
		int db2[2] = {};
		if (fndAddS(db2, resultKeys, the24BitsKey, key1, s1, s2, s4, s5, c1s_2, c2s_2))
			printf("The number of GOST encriptions used in the second stage is: %d\n", NEcr);
		else
		{
			printf("The second stage is failed.\n");
			continue;
		}
		if (!resultKeys.size())
		{
			printf("Failure: not remain keys.\n");
			continue;
		}

		//The third stage.
		int db3[2] = {};
		if (createThirdPCSet(db3, key1, s2, c1s_3, c2s_3, ps))
		{
			printf("The number of GOST encriptions used in the third stage is: %d\n", NEcr);

			vector<uint32_t> resultEntireKeys(256 * resultKeys.size());
			vector<uint64_t> res_s1(256 * resultKeys.size());
			vector<uint64_t> res_s2(256 * resultKeys.size());
			vector<uint64_t> res_s4(256 * resultKeys.size());
			vector<uint64_t> res_s5(256 * resultKeys.size());
			for (int i = 0; i < resultKeys.size(); i++)
			{
				for (uint32_t j = 0; j < 256; j++)
				{
					resultEntireKeys[256 * i + j] = (resultKeys[i] | (j << 24));
					res_s1[256 * i + j] = s1[i];
					res_s2[256 * i + j] = s2[i];
					res_s4[256 * i + j] = s4[i];
					res_s5[256 * i + j] = s5[i];
				}
			}

			int Nk1AfterS7 = 0;
			int Nk1AllSboxes = 0;
			int lastEl[7] = {};
			vector<uint64_t> res_s0(0);
			vector<uint64_t> res_s3(0);
			vector<uint64_t> res_s6(0);
			vector<uint64_t> res_s7(0);
			uint64_t fullK2A[2] = { 0, 0 };
			for (int kind = 0; kind < resultEntireKeys.size(); kind++)
			{
				uint32_t ke = resultEntireKeys[kind];
				int inValsOutDiffs7[NGOODPAIRS][4] = {};

				//Check consistency

				buildInValOutDiffs(7, c1s_3, c2s_3, ke, inValsOutDiffs7, 0, 0, 0);
				int inds7[17] = {};
				if (isWrongKey(inValsOutDiffs7, inds7, 0))
				{
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				int rSboxTemp[8][16] = {};
				buildSbox(rSboxTemp, inValsOutDiffs7, inds7, 7, 0);
				Nk1AfterS7++;
				if (ke == theEntireKey)
				{
					for (int j = 0; j < 16; j++)
						if (rSboxTemp[7][j] ^ s_block[7][j] ^ rSboxTemp[7][0] ^ s_block[7][0])
						{
							printf("The third stage failed.\n");
							continue;
						}
					SUC++;
				}

				//The last stage start here.
				uint64_t c1s_7[2 * NGOODPAIRS] = {};
				uint64_t c2s_7[2 * NGOODPAIRS] = {};
				uint64_t c1s_6[2 * NGOODPAIRS] = {};
				uint64_t c2s_6[2 * NGOODPAIRS] = {};
				uint64_t c1s_4[2 * NGOODPAIRS] = {};
				uint64_t c2s_4[2 * NGOODPAIRS] = {};
				for (int i = 0; i < NGOODPAIRS; i++)
				{
					c1s_7[i] = c1s_3[i];
					c1s_7[i + NGOODPAIRS] = c1s_2[i];
					c1s_4[i] = c1s_2[i];
					c1s_4[i + NGOODPAIRS] = c1s_1[i];
					c1s_6[i] = c1s_3[i];
					c1s_6[i + NGOODPAIRS] = c1s_1[i];
					c2s_7[i] = c2s_3[i];
					c2s_7[i + NGOODPAIRS] = c2s_2[i];
					c2s_4[i] = c2s_2[i];
					c2s_4[i + NGOODPAIRS] = c2s_1[i];
					c2s_6[i] = c2s_3[i];
					c2s_6[i + NGOODPAIRS] = c2s_1[i];
				}

				int inValsOutDiffs0[2 * NGOODPAIRS][4] = {};
				int inValsOutDiffs3[2 * NGOODPAIRS][4] = {};
				int inValsOutDiffs6[2 * NGOODPAIRS][4] = {};
				buildInValOutDiffs(0, c1s_7, c2s_7, ke, inValsOutDiffs0, 1, 0, 0);
				int inds0[17] = {};
				if (isWrongKey(inValsOutDiffs0, inds0, 1))
				{
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				buildInValOutDiffs(3, c1s_4, c2s_4, ke, inValsOutDiffs3, 1, 0, 0);
				int inds3[17] = {};
				if (isWrongKey(inValsOutDiffs3, inds3, 1))
				{
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				buildInValOutDiffs(6, c1s_6, c2s_6, ke, inValsOutDiffs6, 1, 0, 0);
				int inds6[17] = {};
				if (isWrongKey(inValsOutDiffs6, inds6, 1))
				{
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				buildSbox(rSboxTemp, inValsOutDiffs0, inds0, 0, 1);
				buildSbox(rSboxTemp, inValsOutDiffs3, inds3, 3, 1);
				buildSbox(rSboxTemp, inValsOutDiffs6, inds6, 6, 1);
				Nk1AllSboxes++;

				short sb[8][16] = {};
				for (int v = 0; v < 16; v++)
				{
					sb[0][v] = rSboxTemp[0][v];
					sb[1][v] = ((res_s1[kind] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[2][v] = ((res_s2[kind] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[3][v] = rSboxTemp[3][v];
					sb[4][v] = ((res_s4[kind] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[5][v] = ((res_s5[kind] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[6][v] = rSboxTemp[6][v];
					sb[7][v] = rSboxTemp[7][v];
				}

				bool LSB8_A[256] = {};
				bool k2_8[256] = {};

				if (!lastElimination07(c1s_1, c2s_1, sb, ke, LSB8_A, k2_8, DB))
				{
					lastEl[0]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				int nk2 = 0;
				int nlsbA = 0;
				for (int ind8 = 0; ind8 < 256; ind8++)
				{
					if (LSB8_A[ind8])
						nlsbA++;
					if (k2_8[ind8])
						nk2++;
				}
				vector <uint32_t> k2Temp1(nk2);
				int keind = 0;
				for (int i = 0; i < 256; i++)
					if (k2_8[i])
						k2Temp1[keind++] = i;
				vector <uint32_t> ATemp1(nlsbA);
				int Aind = 0;
				for (int i = 0; i < 256; i++)
					if (LSB8_A[i])
						ATemp1[Aind++] = i;

				vector<bool> lsbA_12(nlsbA * 16);
				vector<bool> k2_12(nk2 * 16);
				if (!lastElimination811(c1s_1, c2s_1, sb, ke, lsbA_12, k2_12, DB, ATemp1, k2Temp1))
				{
					lastEl[1]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}
				nk2 = 0;
				nlsbA = 0;
				for (int ind12 = 0; ind12 < lsbA_12.size(); ind12++)
					if (lsbA_12[ind12])
						nlsbA++;
				for (int ind12 = 0; ind12 < k2_12.size(); ind12++)
					if (k2_12[ind12])
						nk2++;
				vector <uint32_t> k2Temp2(nk2);
				keind = 0;
				for (int i = 0; i < k2_12.size(); i++)
					if (k2_12[i])
						k2Temp2[keind++] = (((i % 16) << 8) | (k2Temp1[i / 16]));
				vector <uint32_t> ATemp2(nlsbA);
				Aind = 0;
				for (int i = 0; i < lsbA_12.size(); i++)
					if (lsbA_12[i])
						ATemp2[Aind++] = (((i % 16) << 8) | (ATemp1[i / 16]));

				vector<bool> lsbA_20(nlsbA * 256);
				vector<bool> k2_20(nk2 * 256);
				if (!lastElimination1219(c1s_3, c2s_3, sb, ke, lsbA_20, k2_20, db3, ATemp2, k2Temp2))
				{
					lastEl[2]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}
				nk2 = 0;
				nlsbA = 0;
				for (int ind20 = 0; ind20 < lsbA_20.size(); ind20++)
					if (lsbA_20[ind20])
						nlsbA++;
				for (int ind20 = 0; ind20 < k2_20.size(); ind20++)
					if (k2_20[ind20])
						nk2++;
				k2Temp1.resize(nk2);
				keind = 0;
				for (int i = 0; i < k2_20.size(); i++)
					if (k2_20[i])
						k2Temp1[keind++] = (((i % 256) << 12) | (k2Temp2[i / 256]));
				ATemp1.resize(nlsbA);
				Aind = 0;
				for (int i = 0; i < lsbA_20.size(); i++)
					if (lsbA_20[i])
						ATemp1[Aind++] = (((i % 256) << 12) | (ATemp2[i / 256]));

				vector<bool> lsbA_24(nlsbA * 16);
				vector<bool> k2_24(nk2 * 16);
				if (!lastElimination2023(c1s_3, c2s_3, sb, ke, lsbA_24, k2_24, db3, ATemp1, k2Temp1))
				{
					lastEl[3]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}
				nk2 = 0;
				nlsbA = 0;
				for (int ind24 = 0; ind24 < lsbA_24.size(); ind24++)
					if (lsbA_24[ind24])
						nlsbA++;
				for (int ind24 = 0; ind24 < k2_24.size(); ind24++)
					if (k2_24[ind24])
						nk2++;
				k2Temp2.resize(nk2);
				keind = 0;
				for (int i = 0; i < k2_24.size(); i++)
					if (k2_24[i])
						k2Temp2[keind++] = (((i % 16) << 20) | (k2Temp1[i / 16]));
				ATemp2.resize(nlsbA);
				Aind = 0;
				for (int i = 0; i < lsbA_24.size(); i++)
					if (lsbA_24[i])
						ATemp2[Aind++] = (((i % 16) << 20) | (ATemp1[i / 16]));

				vector<bool> lsbA_32(nlsbA * 256);
				vector<bool> k2_32(nk2 * 256);
				if (!lastElimination2431(c1s_2, c2s_2, sb, ke, lsbA_32, k2_32, db2, ATemp2, k2Temp2))
				{
					lastEl[4]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}
				nk2 = 0;
				nlsbA = 0;
				for (int ind32 = 0; ind32 < lsbA_32.size(); ind32++)
					if (lsbA_32[ind32])
						nlsbA++;
				for (int ind32 = 0; ind32 < k2_32.size(); ind32++)
					if (k2_32[ind32])
						nk2++;
				k2Temp1.resize(nk2);
				keind = 0;
				for (int i = 0; i < k2_32.size(); i++)
					if (k2_32[i])
						k2Temp1[keind++] = (((i % 256) << 24) | (k2Temp2[i / 256]));
				ATemp1.resize(nlsbA);
				Aind = 0;
				for (int i = 0; i < lsbA_32.size(); i++)
					if (lsbA_32[i])
						ATemp1[Aind++] = (((i % 256) << 24) | (ATemp2[i / 256]));
				vector<uint64_t> full_k2_A(0);
				vector<uint32_t> k3Candidates(0);
				if (!lastEliminationR30(c1s_1, c2s_1, c1s_2, c2s_2, c1s_3, c2s_3, sb, ke, DB, db2, db3, ATemp1, k2Temp1, full_k2_A, k3Candidates, theEntireKey, key1[2]))
				{
					lastEl[5]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				uint64_t RBits = 0xffffffff;
				uint64_t LBits = 0xffffffff00000000;

				for (int k2aind = 0; k2aind < full_k2_A.size(); k2aind++)
				{
					vector<uint32_t> k3CandidatesT(k3Candidates.size());
					for (int k3ind = 0; k3ind < k3Candidates.size(); k3ind++)
						k3CandidatesT[k3ind] = k3Candidates[k3ind];
					for (int k3ind = 0; k3ind < k3CandidatesT.size(); k3ind++)
					{
						for (int cind = 0; cind < NGOODPAIRS; cind++)
						{
							uint64_t p1_1 = round(ps[cind], ke, sb);
							p1_1 ^= (full_k2_A[k2aind] & RBits);
							uint64_t p1_2 = round(p1_1, (full_k2_A[k2aind] & LBits) >> 32, sb);
							p1_2 ^= (full_k2_A[k2aind] & RBits);
							uint64_t p1_3 = round(p1_2, k3CandidatesT[k3ind], sb);

							uint64_t p2_1 = round(ps[cind], ke ^ (256 << db3[0]), sb);
							p2_1 ^= (full_k2_A[k2aind] & RBits);
							uint64_t p2_2 = round(p2_1, ((full_k2_A[k2aind] & LBits) >> 32) ^ ((uint32_t)1 << (19 + db3[1])), sb);
							p2_2 ^= (full_k2_A[k2aind] & RBits);
							uint64_t p2_3 = round(p2_2, (k3CandidatesT[k3ind]) ^ (256 << db3[0]), sb);


							if (p1_3^p2_3)
							{
								k3CandidatesT.erase(k3CandidatesT.begin() + (k3ind--));
								break;
							}
						}
					}
					if (!k3CandidatesT.size())
						full_k2_A.erase(full_k2_A.begin() + (k2aind--));
				}
				if (!full_k2_A.size())
				{
					lastEl[6]++;
					resultEntireKeys.erase(resultEntireKeys.begin() + kind);
					res_s1.erase(res_s1.begin() + kind);
					res_s2.erase(res_s2.begin() + kind);
					res_s4.erase(res_s4.begin() + kind);
					res_s5.erase(res_s5.begin() + (kind--));
					continue;
				}

				uint64_t stemp = 0;
				for (int val = 0; val < 16; val++)
					stemp |= ((uint64_t)rSboxTemp[0][val] << (4 * val));
				res_s0.push_back(stemp);
				stemp = 0;
				for (int val = 0; val < 16; val++)
					stemp |= ((uint64_t)rSboxTemp[3][val] << (4 * val));
				res_s3.push_back(stemp);
				stemp = 0;
				for (int val = 0; val < 16; val++)
					stemp |= ((uint64_t)rSboxTemp[6][val] << (4 * val));
				res_s6.push_back(stemp);
				stemp = 0;
				for (int val = 0; val < 16; val++)
					stemp |= ((uint64_t)rSboxTemp[7][val] << (4 * val));
				res_s7.push_back(stemp);
				if (ke == theEntireKey)
				{
					fullK2A[0] = full_k2_A[0];
					fullK2A[1] = full_k2_A[1];
				}
			}

			printf("The number of remain keys after the third stage is: %d\n", Nk1AfterS7);
			printf("The number of remain keys after all the S-boxes are found is: %d\n", Nk1AllSboxes);

			bool fullSboxes = true;
			if (!resultEntireKeys.size())
			{
				fullSboxes = 0;
				printf("The third stage is failed.\in");
				continue;
			}
			for (int ki = 0; ki < resultEntireKeys.size(); ki++)
			{
				short sb[8][16] = {};
				for (int v = 0; v < 16; v++)
				{
					sb[0][v] = ((res_s0[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[1][v] = ((res_s1[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[2][v] = ((res_s2[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[3][v] = ((res_s3[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[4][v] = ((res_s4[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[5][v] = ((res_s5[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[6][v] = ((res_s6[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
					sb[7][v] = ((res_s7[ki] & (((uint64_t)0xf) << (4 * v))) >> (4 * v));
				}
				uint32_t sbxor0 = fullK2A[0];
				uint32_t sbxor1 = fullK2A[1];
				sbxor0 = (sbxor0 >> 11) | (sbxor0 << 22);
				sbxor1 = (sbxor1 >> 11) | (sbxor1 << 22);
				printf("The S-box for the subkey K1 option %lx is one of the following two options:\n", resultEntireKeys[ki]);
				for (int sbx = 0; sbx < 8; sbx++)
				{
					for (int sv = 0; sv < 16; sv++)
						printf("%x\t", (sb[sbx][sv] ^ ((sbxor0&(0xf << 4 * sbx)) >> 4 * sbx)));
					printf("\n");
				}
				printf("Or:\n");
				for (int sbx = 0; sbx < 8; sbx++)
				{
					for (int sv = 0; sv < 16; sv++)
						printf("%x\t", (sb[sbx][sv] ^ ((sbxor1&(0xf << 4 * sbx)) >> 4 * sbx)));
					printf("\n");
				}
				printf("The 2 options for the key K2 are: %lx\t%lx\n",
					(fullK2A[0] >> 32), (fullK2A[1] >> 32));
				if (resultEntireKeys[ki] == theEntireKey)
				{
					for (int s = 0; s < 8; s++)
					{
						for (int v = 1; v < 16; v++)
						{
							if (sb[s][v] ^ s_block[s][v] ^ sb[s][0] ^ s_block[s][0])
								fullSboxes = false;
						}
					}
				}
			}

			if (fullSboxes)
			{
				if ((((fullK2A[0] & 0x780) >> 7) == s_block[7][0]) ||
					(((fullK2A[1] & 0x780) >> 7) == s_block[7][0]))
				{
					SUCSUC++;
					printf("%d out of %d experiments was succeeded.\n",
						SUCSUC, expInd + 1);
				}
				else
				{
					SUCSUCOI++;
					cout << "The recovered S-boxes are wrong!" << endl;
				}

				printf("The number of remain K1 is: %d\n", resultEntireKeys.size());
			}
		}
		else
			printf("The third stage is failed.\n");
	}
	printf("Time = %d\n", time(0));
	printf("The success rate of the 1st, 2nd, 3rd steps are: %d\t%d\t%d\n",
		SUC1, SUC2, SUC);
	printf("The success rate of the entire attack is: %d\n", SUCSUC);
}

uint64_t MDGost(vector<uint64_t> M, uint64_t lenM, uint64_t IV)
{
	//Padding steps
	if (lenM % 256)
		M[lenM >> 6] |= ((uint64_t)1 << (lenM % 64));	//Add 1 after
	else
		M.push_back(1);

	//Add zeros
	uint64_t numOfZeros = ((BLOCKSIZE - ((lenM + 1 + MAXMESSAGELEN) % BLOCKSIZE)) % BLOCKSIZE);
	if (numOfZeros > 63)
	{
		uint64_t n = numOfZeros >> 6;
		while (n)
		{
			M.push_back(0);
			n--;
		}
	}
	//Add the message length
	M.push_back(lenM);

	if (M.size() % 4)
		printf("problem!\n");

	//perform the MD construction

	int numOfIter = M.size() >> 2;

	uint64_t block = IV;
	for (int i = 0; i < numOfIter; i++)
	{
		uint32_t key[8] = {};
		for (int k = 0; k < 4; k++)
		{
			uint64_t lsb = M[4 * i + k] & 0xffffffff;
			uint64_t msb = M[4 * i + k] >> 32;
			key[2 * k] = lsb;
			key[2 * k + 1] = msb;
		}
		block = encrypt(block, key, s_block, 32);
	}

	return block;
}

int chooseDiffBit(double sump[4])
{
	int ddt[16][16] = {};
	calcDDT(s_block[7], ddt);

	int diffBit = 0;

	sump[0] = (((double)ddt[8][1] / 16.))
		*(((double)ddt[8][1] / 16.))*0.5;
	sump[0] += (((double)ddt[8][3] / 16.))
		*(((double)ddt[8][3] / 16.))*0.25;
	sump[0] += (((double)ddt[8][7] / 16.))
		*(((double)ddt[8][7] / 16.))*0.125;
	sump[0] += (((double)ddt[8][15] / 16.))
		*(((double)ddt[8][15] / 16.))*0.0625;

	sump[1] = (((double)ddt[8][2] / 16.))
		*(((double)ddt[8][2] / 16.))*0.5;
	sump[1] += (((double)ddt[8][6] / 16.))
		*(((double)ddt[8][6] / 16.))*0.25;
	sump[1] += (((double)ddt[8][14] / 16.))
		*(((double)ddt[8][14] / 16.))*0.125;

	sump[2] = (((double)ddt[8][4] / 16.))
		*(((double)ddt[8][4] / 16.))*0.5;
	sump[2] += (((double)ddt[8][12] / 16.))
		*(((double)ddt[8][12] / 16.))*0.25;

	sump[3] = (((double)ddt[8][8] / 16.))
		*(((double)ddt[8][8] / 16.))*0.5;

	double maxp = sump[0];
	for (int i = 0; i < 4; i++)
	{
		if (sump[i] > maxp)
		{

			maxp = sump[i];
			diffBit = i;
		}
	}

	return diffBit;
}

int fndM(uint64_t inpt, uint32_t M1[8], uint32_t M2[8])
{
	double sump[4] = {};
	int db = chooseDiffBit(sump);
	bool notFound = 1;
	int numEx = 0;
	while (notFound)
	{
		numEx++;
		if (numEx == 1000000)
			break;
		uint32_t key1[8] = {};
		for (int i = 0; i < 6; i++)
			key1[i] = ((uint64_t)rand()) | ((uint64_t)(rand()) << 15) | ((uint64_t)(rand() % 4) << 30);
		uint32_t key2[8] = {};
		for (int i = 0; i < 8; i++)
			key2[i] = key1[i];
		key2[0] ^= ((uint32_t)1 << 31);
		key2[1] ^= ((uint32_t)1 << (7 + db));
		key2[2] ^= ((uint32_t)1 << 31);

		uint64_t X1 = encrypt(inpt, key1, s_block, 3);
		uint64_t X2 = encrypt(inpt, key2, s_block, 3);

		if (X1 == X2)
		{
			notFound = 0;
			X1 = encrypt(inpt, key1, s_block, 6);
			X2 = encrypt(inpt, key2, s_block, 6);
			X1 = (X1 << 32) | (X1 >> 32);
			X2 = (X2 << 32) | (X2 >> 32);
			uint64_t xr = X1 & 0xffffffff;
			uint64_t xl = X1 >> 32;
			uint64_t pr = inpt & 0xffffffff;
			uint64_t pl = inpt >> 32;
			uint64_t sInvR = 0;
			uint32_t bfAfXor = pr ^ xr;
			uint32_t mask = bfAfXor >> 11;
			bfAfXor = (bfAfXor << 21) | mask;
			for (int i = 0; i < 8; i++)
				sInvR |= (((uint64_t)sInverse[i][(bfAfXor&(0xf << 4 * i)) >> 4 * i]) << 4 * i);
			uint64_t sInvL = 0;
			bfAfXor = pl ^ xl;
			mask = bfAfXor >> 11;
			bfAfXor = (bfAfXor << 21) | mask;
			for (int i = 0; i < 8; i++)
				sInvL |= (((uint64_t)sInverse[i][(bfAfXor&(0xf << 4 * i)) >> 4 * i]) << 4 * i);
			uint64_t twoto32 = 0x100000000;
			key1[6] = sInvL - xr;
			key1[7] = sInvR - pl;

			//key1[6] = ((twoto32 + sInvL - xr) % twoto32);
			//key1[7] = ((twoto32 + sInvR - pl) % twoto32);
			key2[6] = key1[6];
			key2[7] = key1[7];

			X1 = encrypt(inpt, key1, s_block, 32);
			X2 = encrypt(inpt, key2, s_block, 32);
			int d = 1;
			while (X1^X2)
			{
				d++;
				if (d == 1000000)
					break;
				key1[5] ^= d;
				key2[5] = key1[5];

				X1 = encrypt(inpt, key1, s_block, 6);
				X2 = encrypt(inpt, key2, s_block, 6);
				X1 = (X1 << 32) | (X1 >> 32);
				X2 = (X2 << 32) | (X2 >> 32);
				uint64_t xr = X1 & 0xffffffff;
				uint64_t xl = X1 >> 32;
				uint64_t pr = inpt & 0xffffffff;
				uint64_t pl = inpt >> 32;
				uint64_t sInvR = 0;
				uint32_t bfAfXor = pr ^ xr;
				uint32_t mask = bfAfXor >> 11;
				bfAfXor = (bfAfXor << 21) | mask;
				for (int i = 0; i < 8; i++)
					sInvR |= (((uint64_t)sInverse[i][(bfAfXor&(0xf << 4 * i)) >> 4 * i]) << 4 * i);
				uint64_t sInvL = 0;
				bfAfXor = pl ^ xl;
				mask = bfAfXor >> 11;
				bfAfXor = (bfAfXor << 21) | mask;
				for (int i = 0; i < 8; i++)
					sInvL |= (((uint64_t)sInverse[i][(bfAfXor&(0xf << 4 * i)) >> 4 * i]) << 4 * i);
				uint64_t twoto32 = 0x100000000;
				key1[6] = sInvL - xr;
				key1[7] = sInvR - pl;
				key2[6] = key1[6];
				key2[7] = key1[7];

				X1 = encrypt(inpt, key1, s_block, 32);
				X2 = encrypt(inpt, key2, s_block, 32);
				if (X1 == X2)
				{
					notFound = 0;
					for (int a = 0; a < 8; a++)
					{
						M1[a] = key1[a];
						M2[a] = key2[a];
					}
				}
			}
			numEx += d;
		}
	}
	if (notFound)
		return -1;
	return numEx;
}

void CollisionInMDGost()
{
	int sumEx = 0;
	int suc = 0;
	for (int expInd = 0; expInd < 100; expInd++)
	{
		printf("Experiment num. %d\n", expInd);
		for (int s = 0; s < 8; s++)
		{
			int tempArr[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
			randomShuffle(tempArr, 16);
			for (int val = 0; val < 16; val++)
				s_block[s][val] = tempArr[val];
		}
		initSInv();

		uint64_t iv = 0;
		for (int i = 0; i < 4; i++)
			iv |= (((uint64_t)rand()) << (15 * i));
		iv |= ((((uint64_t)rand()) % 16) << 60);

		uint32_t M1[8] = {};
		uint32_t M2[8] = {};
		int nEncryptions = fndM(iv, M1, M2);
		if (nEncryptions > -1)
			sumEx += nEncryptions;
		else
		{
			printf("Failure");
			continue;
		}
		vector<uint64_t> m1(0);
		m1.resize(4);
		for (int i = 0; i < 4; i++)
			m1[i] = (((uint64_t)M1[2 * i + 1]) << 32) | (M1[2 * i]);
		vector<uint64_t> m2(0);
		m2.resize(4);
		for (int i = 0; i < 4; i++)
			m2[i] = (((uint64_t)M2[2 * i + 1]) << 32) | (M2[2 * i]);
		uint64_t output1 = MDGost(m1, 256, iv);
		uint64_t output2 = MDGost(m2, 256, iv);
		if (output1 == output2)
		{
			suc++;
			printf("Let IV be %llx.\n", iv);
			printf("The GOST-based MD of the message\n");
			for (int i = 0; i < 4; i++)
				printf("%llx\t", m1[i]);
			printf("\n");
			printf("is %llx\n", output1);
			printf("The GOST-based MD of the message\n");
			for (int i = 0; i < 4; i++)
				printf("%llx\t", m2[i]);
			printf("\n");
			printf("is %llx\n", output2);
			printf("The collision is found by %d GOST encryptions.\n", nEncryptions);
		}
		else
			printf("Failure");
	}
	double avg = ((double)sumEx) / ((double)suc);
	printf("The success rate is: %d\n", suc);
	printf("The avg. of GOST encryptions is: %f\n", avg);
}

void example()
{
	//Example on the bank S-boxes
	for (int s = 0; s < 8; s++)
		for (int v = 0; v < 16; v++)
			s_block[s][v] = bankS[s][v];
	initSInv();
	int sumEx = 0;
	uint32_t M1[8] = {};
	uint32_t M2[8] = {};
	int nEncryptions = fndM(0, M1, M2);
	if (nEncryptions > -1)
	{
		sumEx += nEncryptions;
		vector<uint64_t> m1(0);
		m1.resize(4);
		for (int i = 0; i < 4; i++)
			m1[i] = (((uint64_t)M1[2 * i + 1]) << 32) | (M1[2 * i]);
		vector<uint64_t> m2(0);
		m2.resize(4);
		for (int i = 0; i < 4; i++)
			m2[i] = (((uint64_t)M2[2 * i + 1]) << 32) | (M2[2 * i]);
		uint64_t output1 = MDGost(m1, 256, 0);
		uint64_t output2 = MDGost(m2, 256, 0);

		printf("Example1: Collision for the Banking industry S-boxes and IV=0:\n");
		printf("M1: \n");
		for (int i = 0; i < 4; i++)
			printf("%llx\t", m1[i]);
		printf("\n");
		printf("Output: %llx\n", output1);
		printf("M2: \n");
		for (int i = 0; i < 4; i++)
			printf("%llx\t", m2[i]);
		printf("\n");
		printf("Output: %llx\n", output2);
	}

	else
		printf("Failure");

	//Example on the gost2 S-boxes
	for (int s = 0; s < 8; s++)
		for (int v = 0; v < 16; v++)
			s_block[s][v] = gost2S[s][v];
	initSInv();
	sumEx = 0;
	M1[8] = {};
	M2[8] = {};
	nEncryptions = fndM(0, M1, M2);
	if (nEncryptions > -1)
	{
		sumEx += nEncryptions;
		vector<uint64_t> m1(0);
		m1.resize(4);
		for (int i = 0; i < 4; i++)
			m1[i] = (((uint64_t)M1[2 * i + 1]) << 32) | (M1[2 * i]);
		vector<uint64_t> m2(0);
		m2.resize(4);
		for (int i = 0; i < 4; i++)
			m2[i] = (((uint64_t)M2[2 * i + 1]) << 32) | (M2[2 * i]);
		uint64_t output1 = MDGost(m1, 256, 0);
		uint64_t output2 = MDGost(m2, 256, 0);

		printf("Example2: Collision for the GOST2 S-boxes and IV=0:\n");
		printf("M1: \n");
		for (int i = 0; i < 4; i++)
			printf("%llx\t", m1[i]);
		printf("\n");
		printf("Output: %llx\n", output1);
		printf("M2: \n");
		for (int i = 0; i < 4; i++)
			printf("%llx\t", m2[i]);
		printf("\n");
		printf("Output: %llx\n", output2);
	}

	else
		printf("Failure");
}

int main()
{
	time_t ttt;
	srand((unsigned)time(&ttt));
	uint16_t pp = 0x7fff;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0, pp);
	initSInv();;
	printf("The experiment of Section 4 begins now.");
	keySboxRecoveryAttackZD();
	printf("The experiment of Section 5.1 begins now.");
	CollisionInMDGost();
	example();

	return 0;
}
