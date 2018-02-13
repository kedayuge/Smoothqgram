#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <tuple>
#include <omp.h>
#include <numeric>
#include <math.h>
#include <map>
#include <set>
#include <climits>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <iterator>
#include <thread>
#include "xxHash/xxhash.c"
#include "edit.h"
#include <time.h>
#include <sys/time.h>
using namespace std;

typedef tuple <int, int> intpair;

struct signature
{
	int location, id;
};

struct signature2
{
	int location;
	uint64_t value;
};

struct match
{
	int loc1, loc2;
};

struct Data {
	string str; //initialize in readdata()
	const char * cstr;
	int len;
	int index; // indices of original data
	vector<string> ebdqgrams;

	int fwrv; // fwd=0 rev=1;
	vector<signature2> sig;
};

vector<Data> oridata; //input data
unordered_set<string> freqqgrams;// frequenct smooth q-gram to be filtered 
string filename; //input file name

//for running time of verification
double vtime = 0.0;

// for DEBUG
int DEBUG = 0;
int id1 = 0, id2 = 0, fwd = 0;

unordered_map<char, char> trans({ { 'A', 'T' }, { 'T', 'A' }, { 'G', 'C' }, { 'C', 'G' } });
unordered_map<char, int> dict({ { 'A', 0 }, { 'C', 1 }, { 'G', 2 }, { 'T', 3 }, { 'N', 4 } });


//Parameters to be tuned

int max_data = 2000000; // maximum number of input string
int min_len = 200; // minimum length of input to consider

int q = 14; //len of qgram
double samplerate = 1.5; // rate for random sampling, len of smooth q-gram = samplerate * q

double hashpercent = 0.15;// signature selection rate
double freqthres = 0.00003; //frequency filtering threshold

int windowsize = 500; //targeting overlap length, minimum overlap length to output

//parameters for verification
double threshold = 3; // threshold for #matched sigs
double shiftpercent = 0.2; // error rate of SMRT
int extendstep = 1000; // maximum step to extend in verification




