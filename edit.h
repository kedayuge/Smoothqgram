#pragma once
// A C++ program to find k'th largest element in a stream

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

using namespace std;

#define K  (4)
typedef int64_t int64;
typedef int32_t int32;

int slide(const char *x, const char *y);

int slide32(const char *x, const char *y);

int edit_dp(const char *x, const int x_len, const  char *y, const int y_len, int k);


