#include "pch.h"


const int LATTICE_SIZE = 8;
const int NUM_SITES = LATTICE_SIZE * LATTICE_SIZE;
const int NUM_CONFIGS = NUM_SITES / 32;


// ls = lattice size

// returns a mod b, where b is positive
int mod(int a, int b) {
	int ret = a % b;
	if (ret < 0) {
		ret += b;
		return mod(ret, b);
	}
	return ret;
}


int int32rand() {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, 65535);
	// juxtapose two random short to get a random unsigned int 
	return (dis(gen) << 16) | dis(gen);
	// the final displayed number depends on how the OS explains the string of 101000111... bits
	// the string of bits is fixed and sufficient for our purpose
}


// get the bit at position (i,j) of the lattice represented by A[2]
int get_bit(int i, int j, int* A) {
	int k = i * LATTICE_SIZE + j;
	int index = k / 32;
	int bit_pos = k % 32;
	int n = *(A + index);

	//return (n&(1<< bit_pos)) != 0 ;
	return (n >> bit_pos) & 1;
}

// int** a is a pointer to int A[2]
// usage: 
// int* p = A;
//set_bit(i, j, s, &p)
// or direct set_bit(i, j, s, &A)
void set_bit(int i, int j, int spin_value, int** a) {
	int k = i * LATTICE_SIZE + j;
	int index = k / 32;
	int bit_pos = k % 32;

	if (spin_value == 1) {
		*a[index] = *a[index] | (1 << bit_pos);
	}
	else {
		*a[index] = *a[index] & ~(1 << bit_pos);
	}
}


// nearest neighbour check
int nn_check(int i1, int j1, int i2, int j2) {
	int x_diff = mod(i1 - j1, LATTICE_SIZE);
	int y_diff = mod(i2 - j2, LATTICE_SIZE);
	if ((abs(x_diff) == 1 && (y_diff == 0)) ||
		(abs(y_diff) == 1 && (x_diff == 0))) return 1;
	else return 0;
}

// calculate energy of a spin pair
/*
int pair_energy(int i1, int j1, int i2, int j2, int* A){
	if (! nn_check(i1, j1, i2, j2)) printf("%s\n", "Error! Non-pair encountered.");
	else {
		if (get_bit(i1, j1, A) == get_bit(i2, j2, A)) return 2;
		else return -2;
	}
}*/

int pair_energy(int s1, int s2) {
	return (s1 ^ s2) ? -2 : 2;
}


// calculate the energy cost of a spin flip at position (i,j)
int dE(int i, int j, int* A) {
	return pair_energy(get_bit(i, j, A), !get_bit(i, mod(j - 1, LATTICE_SIZE), A)) + \
		pair_energy(get_bit(i, j, A), !get_bit(i, mod(j + 1, LATTICE_SIZE), A)) + \
		pair_energy(get_bit(i, j, A), !get_bit(mod(i - 1, LATTICE_SIZE), j, A)) + \
		pair_energy(get_bit(i, j, A), !get_bit(mod(i + 1, LATTICE_SIZE), j, A))
		;
}


//calculate the magnetization for a given configuration A
//which is equal to the count of all bits that are set to 1
int bit_count(int n) {
	int count = 0;
	while (n != 0) {
		if (n & 1) count++;
		n = (n >> 1);
	}
	return count;
}

int magnetization(int* A) {
	int mag = 0;
	for (int i = 0; i < NUM_CONFIGS; ++i) {
		mag += bit_count(A[i]);
	}
	return 2 * mag - NUM_SITES;
}


int initialization()
{
	// initialization
	int A[2] = { int32rand(), int32rand() };

	// test
	printf("%d\n", A[0]);
	printf("%d\n", A[1]);
	for (int i = 0; i < LATTICE_SIZE; ++i) {
		for (int j = 0; j < LATTICE_SIZE; ++j)
			printf("%d ", get_bit(i, j, A));
	}
	return 0;
}



int main() {
	unsigned int i = -7;
	printf("%d", (i));
	return 0;
}