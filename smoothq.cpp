#include "smoothq.h"


inline string qgrams(int i, int j)
{
	return oridata[i].str.substr(j, q);
}

double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}

inline bool sortintpair(const intpair& a, const intpair& b)
{
	// smallest comes first
	return get<0>(a) < get<0>(b);
}

inline bool sortsig2(const signature2& a, const signature2& b)
{
	// smallest comes first
	return a.value < b.value;
}

inline bool sortmatch(const match& a, const match& b)
{
	// smallest comes first
	return a.loc2 - a.loc1 < b.loc2 - b.loc1;
}

inline bool sortmatch2(const match& a, const match& b)
{
	// smallest comes first
	return a.loc1 < b.loc1;
}

template <class InIt, class OutIt, class size_type>
void window(InIt b, InIt e, OutIt d, size_type len) {
	for (auto s = std::next(b, len); s != e; ++s, ++b, ++d)
		*d = { b, s };
}

vector<string> generate_qgram(string const &in, size_t q) {
	string input = in;
	vector<string> words;
	window(input.begin(), input.end(), back_inserter(words), q);
	return words;
}

string revstring(string& input)
{
	string output;
	output.resize(input.size());

	transform(input.rbegin(), input.rend(), output.begin(), [](char c) -> char {
		return trans[c];
	});

	return output;
}

//load data and generate qgrams
void readfasta()
{
	ifstream input(filename);
	string line, str;
	int avglen = 0;
	int maxlen = 0;
	int index = 1;

	while (getline(input, line) && oridata.size() < max_data) {
		if (line.empty())
			continue;
		if (line[0] == '>')
		{
			str.erase(remove(str.begin(), str.end(), 13), str.end());
			if (str != "")
			{
				if (str.size() > min_len)
				{
					// add the forward sequence
					Data newdata;
					newdata.str = str;
					newdata.index = index;
					newdata.len = str.size();
					newdata.fwrv = 0;
					newdata.cstr = newdata.str.c_str();
					oridata.push_back(newdata);

					// add the reverse sequence
					newdata.str = revstring(str);
					newdata.fwrv = 1;
					newdata.cstr = newdata.str.c_str();
					oridata.push_back(newdata);

					avglen += str.size();
					if (str.size() > maxlen)
						maxlen = str.size();
				}
				index++;
			}
			str.clear();
		}
		else
		{
			str += line;
		}
	}

	fprintf(stderr, "Load %d sequences\n", oridata.size());
	fprintf(stderr, "The average length of sequence is: %d\n", avglen / oridata.size());
	fprintf(stderr, "The maximum length of sequence is: %d\n", maxlen);
}

void countqgramfreq(unordered_set<string>& output, double threshold)
{
	unordered_map<string, int> freq;
	for (auto & point : oridata)
		for (auto & qgram : point.ebdqgrams)
			++freq[qgram];

	// record frequent qgrams
	for (auto & x : freq)
		if (x.second >= threshold * freq.size()) //x.second > 0.000005*freq.size() &&
		{
			output.insert(x.first);
		}
}

inline string embedstr(string& input, int **p, vector<int>& hash_lsh)
{
	string output;
	int i = 0, partdigit = 0;

	for (int j = 0; partdigit < hash_lsh.size(); ++j)
	{
		char s = i < input.size() ? input[i] : 'N';
		//only record the bits used in LSH
		while (partdigit < hash_lsh.size() && j == hash_lsh[partdigit])
		{
			output += s;
			partdigit += 1;
		}
		i = i + *(p[dict[s]] + j);
	}

	return output;
}

vector<string> embeddata(int i, int** p, vector<int>& hash_lsh)
{
	vector<string> qgrams = generate_qgram(oridata[i].str, q);
	vector<string> output;
	output.reserve(qgrams.size());

	for (int i = 0; i < qgrams.size(); ++i)
	{
		output.emplace_back(embedstr(qgrams[i], p, hash_lsh));
	}

	return output;
}

vector<int> producelsh(int q)
{
	vector<int> hash_lsh;
	for (int j = 0; j < 2 * q; j++)
		hash_lsh.push_back(j);
	random_shuffle(hash_lsh.begin(), hash_lsh.end());
	hash_lsh.resize(q*samplerate);
	//hash_lsh.push_back(rand() % (2 * q));
	sort(hash_lsh.begin(), hash_lsh.end()); // remove the duplicates

	return hash_lsh;
}

int** produceran(int q)
{
	const int num_char = 5;
	int **p = new int *[num_char];
	int *hash_eb = new int[num_char * (2 * q + 1)];

	for (int t = 0; t < num_char; t++)
	{
		p[t] = &hash_eb[t*(2 * q + 1)];
		for (int d = 0; d < 2 * q + 1; d++)
		{
			hash_eb[t*(2 * q + 1) + d] = rand() % 2;
		}
	}
	return p;
}

void qgramandhash()
{
	int seed = 1702;

	vector<int> hash_lsh = producelsh(q);
	int **p = produceran(q);

	#pragma omp parallel for
	for (int i = 0; i < oridata.size(); i++)
	{
		oridata[i].ebdqgrams = embeddata(i, p, hash_lsh);
	}
	countqgramfreq(freqqgrams, freqthres);

	#pragma omp parallel for
	for (int i = 0; i < oridata.size(); i++)
	{
		// hash qgrams
		vector<signature2> sig2;
		for (int j = 0; j < oridata[i].ebdqgrams.size(); j++)
			if (freqqgrams.find(oridata[i].ebdqgrams[j]) == freqqgrams.end())
			{
				const char *cstr = oridata[i].ebdqgrams[j].c_str();
				sig2.push_back({ j, XXH32(cstr, strlen(cstr), seed) });
			}

		sort(sig2.begin(), sig2.end(), sortsig2);
		if (oridata[i].len> sig2.size())
			sig2.resize(oridata[i].len);
		oridata[i].sig = sig2;
	}

	fprintf(stderr, "finish hash strings \n");
}

unordered_map<string, list<signature>> buildtable()
{
	unordered_map<string, list<signature>> hash;

	for (int i = 0; i < oridata.size(); i++) //id
	{
		int count = 0;
		// for each hashed value
		for (int j = 0; count < oridata[i].len*hashpercent; ++j)
		{
			hash[oridata[i].ebdqgrams[oridata[i].sig[j].location]].push_back({ oridata[i].sig[j].location,i });
			++count;
		}
	}

	fprintf(stderr, "finish build hash table \n");
	return hash;
}

void findoffset(vector<match>& matches, int maxshift, vector<match> & validmatches, int DEBUG)
{
	int maxcount = 0, maxmedian = 0;
	// record count of offset of matches
	vector<intpair> offsets; // offset, start/end(0/1)
	for (auto& x : matches)
	{
		offsets.push_back(make_tuple(x.loc2 - x.loc1 - maxshift, 0));
		offsets.push_back(make_tuple(x.loc2 - x.loc1 + maxshift, 1));
	}
	sort(offsets.begin(), offsets.end(), sortintpair);

	int count = 0;
	for (auto& x : offsets)
	{
		if (get<1>(x) == 0)
		{
			count++;
			if (count > maxcount)
			{
				maxcount = count;
				maxmedian = get<0>(x);
			}
		}
		else
		{
			count--;
		}
	}

	for (auto& x : matches)
	{
		if (abs(x.loc2 - x.loc1 - maxmedian) <= maxshift)
		{
			validmatches.push_back(x);
		}
	}

	if (DEBUG)
	{
		sort(matches.begin(), matches.end(), sortmatch2);
		sort(validmatches.begin(), validmatches.end(), sortmatch2);
		cout << "The matches are:" << endl;
		for (auto x : matches)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
		cout << "The valid matches are:" << endl;
		for (auto x : validmatches)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
	}
}

int findoverlaparea(vector<match>& matches, vector<match> & validmatches, int DEBUG)
{
	vector<intpair> locs; // offset, start/end(0/1)
	for (auto& x : matches)
	{
		locs.push_back(make_tuple(x.loc1 - windowsize / 2, 0));
		locs.push_back(make_tuple(x.loc1 + windowsize / 2, 1));
	}
	sort(locs.begin(), locs.end(), sortintpair);

	int maxcount = 0, xloc = 0;
	int count = 0;
	for (auto& x : locs)
	{
		if (get<1>(x) == 0)
		{
			count++;
			if (count > maxcount)
			{
				maxcount = count;
				xloc = get<0>(x);
			}
		}
		else
		{
			count -= 1;
		}
	}

	for (auto& x : matches)
	{
		if (abs(x.loc1 - xloc) <= windowsize / 2)
			validmatches.push_back(x);
	}
	sort(validmatches.begin(), validmatches.end(), sortmatch2);
	if (DEBUG)
	{
		cout << "The valid matches in the window are:" << endl;
		for (auto x : validmatches)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
		cout << "Count:" << maxcount << " with x starting at:" << xloc + windowsize / 2 << endl;
	}

	return maxcount;
}

void extendmatch(vector<match>& matches, vector<match> & extendedmatches, int DEBUG)
{
	vector<match> prematch, backmatch;
	sort(matches.begin(), matches.end(), sortmatch2);
	sort(extendedmatches.begin(), extendedmatches.end(), sortmatch2);
	match xbegin = extendedmatches[0], xend = extendedmatches[extendedmatches.size() - 1];

	for (auto& x : matches)
	{
		if (x.loc1 < xbegin.loc1)
		{
			prematch.push_back(x);
		}
		if (x.loc1 > xend.loc1)
		{
			backmatch.push_back(x);
		}
	}
	sort(prematch.begin(), prematch.end(), sortmatch2);
	reverse(prematch.begin(), prematch.end());
	sort(backmatch.begin(), backmatch.end(), sortmatch2);

	for (int i = 0; i < prematch.size() && (xbegin.loc1 - prematch[i].loc1 < extendstep); i++)
	{
		int offsetshift = abs(prematch[i].loc2 - prematch[i].loc1 - xbegin.loc2 + xbegin.loc1);
		int locshift = max(windowsize, abs(prematch[i].loc1 - xbegin.loc1));
		if (offsetshift <= locshift*shiftpercent)
		{
			xbegin = prematch[i];
			extendedmatches.push_back(xbegin);
		}
	}

	for (int i = 0; i < backmatch.size() && (backmatch[i].loc1 - xend.loc1 < extendstep); i++)
	{
		int offsetshift = abs(backmatch[i].loc2 - backmatch[i].loc1 - xend.loc2 + xend.loc1);
		int locshift = max(windowsize, abs(backmatch[i].loc1 - xend.loc1));
		if (offsetshift <= locshift*shiftpercent)
		{
			xend = backmatch[i];
			extendedmatches.push_back(xend);
		}
	}

	if (DEBUG)
	{
		sort(extendedmatches.begin(), extendedmatches.end(), sortmatch2);
		cout << "The preMatches are:" << endl;
		for (auto x : prematch)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
		cout << "The backmatches are:" << endl;
		for (auto x : backmatch)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
		cout << "The extended valid matches are:" << endl;
		for (auto x : extendedmatches)
		{
			cout << x.loc1 << " " << x.loc2 << " " << x.loc2 - x.loc1 << endl;
		}
	}
}

void addoutput(int i, int j, double maxcount, vector<match>& validmatches)
{
	double maxshift = windowsize * shiftpercent;
	int a1, a2, b1, b2;
	a2 = b2 = 0;
	a1 = oridata[i].len;
	b1 = oridata[j].len;

	for (auto x : validmatches)
	{
		if (x.loc1 < a1)
			a1 = x.loc1;
		if (x.loc1 > a2)
			a2 = x.loc1;

		if (x.loc2 < b1)
			b1 = x.loc2;
		if (x.loc2 > b2)
			b2 = x.loc2;
	}

	if (oridata[j].fwrv)
	{
		int tmp = b2;
		b2 = oridata[j].len - b1;
		b1 = oridata[j].len - tmp;
	}

	int newa1, newa2, newb1, newb2;
	newa1 = max(0, (int)((double)(a1*validmatches.size() - a2) / (double)(validmatches.size() - 1)));
	newa2 = min((int)oridata[i].len, (int)((double)(a2*validmatches.size() - a1) / (double)(validmatches.size() - 1)));
	newb1 = max(0, (int)((double)(b1*validmatches.size() - b2) / (double)(validmatches.size() - 1)));
	newb2 = min((int)oridata[i].len, (int)((double)(b2*validmatches.size() - b1) / (double)(validmatches.size() - 1)));
	/*
	if (newa2 - newa1>100)
	cout << oridata[i].index << " " << oridata[j].index << " " << 0.1 << " " << maxcount << " " << oridata[i].fwrv << " " << newa1 << " " << newa2 << " " << oridata[i].len << " " << oridata[j].fwrv << " " << newb1 << " " << newb2 << " " << oridata[j].len << endl;
	*/
	if (a2 - a1 > 100)
		cout << oridata[i].index << " " << oridata[j].index << " " << maxcount << " " << oridata[i].fwrv << " " << a1 << " " << a2 << " " << oridata[i].len
		<< " " << oridata[j].fwrv << " " << b1 << " " << b2 << " " << oridata[j].len << endl;
}

vector<match> calmatches(int i, int j, int DEBUG)
{
	vector<signature2> x = oridata[i].sig, y = oridata[j].sig;
	vector<match> output;
	int repeat = 1;
	// init counters
	int i1 = 0;
	int i2 = 0;
	// perform merge operation to get the shift and the kmer count

	if (DEBUG)
	{
		cout << "matches in the second stage:" << endl;
	}

	while (true)
	{
		if (i1 >= x.size())
			break;
		if (i2 >= y.size())
			break;
		// get the values in the array
		uint64_t hash1 = x[i1].value;
		uint64_t hash2 = y[i2].value;

		if (DEBUG)
		{
			//cout << hash1 << " " << hash2 << endl;
		}

		if (hash1 < hash2)
			i1++;
		else if (hash2 < hash1)
			i2++;
		else
		{
			//record match
			if (edit_dp(&(oridata[i].cstr[x[i1].location]), q, &(oridata[j].cstr[y[i2].location]), q, 2) != -1)
			{
				output.push_back({ x[i1].location, y[i2].location });
			}
			if (repeat == 0)
				i1++;
			i2++;
		}
	}
	return output;
}

int verify(int i, int j, vector<match>& matches, int DEBUG)
{
	if (DEBUG)
	{
		cout << oridata[i].str << endl << oridata[j].str << endl;
	}
	if (!DEBUG && matches.size() < threshold)
		return 0;

	double maxshift = windowsize*shiftpercent;
	vector<match> validmatches, validmatches2;
	findoffset(matches, maxshift, validmatches, DEBUG);
	if (!DEBUG && validmatches.size() < threshold)
		return 0;
	int count = findoverlaparea(validmatches, validmatches2, DEBUG);

	if (count >= threshold && (abs(validmatches2[0].loc1 - validmatches2[validmatches2.size() - 1].loc1) > 100))
	{
		int off = validmatches2[0].loc1 - validmatches2[0].loc2, overlen;
		if (off > 0)
			overlen = min(oridata[i].len - off, oridata[j].len);
		else
			overlen = min(oridata[j].len + off, oridata[i].len);
		vector<match> newmatches, newvalidmatches, newvalidmatches2;
		newmatches = calmatches(i, j, DEBUG);
		findoffset(newmatches, max(maxshift, maxshift), newvalidmatches, DEBUG);
		if (newvalidmatches.size() >= threshold)
		{
			findoverlaparea(newvalidmatches, newvalidmatches2, DEBUG);
			extendmatch(newmatches, newvalidmatches2, DEBUG);
			addoutput(i, j, count, newvalidmatches2);
		}
	}
}

void join(unordered_map<string, list<signature>>& hash)
{
	for (int i = 0; i < oridata.size(); i++)
	{
		unordered_map<int, vector<match>> candidate;
		//record the candidate
		for (int j = 0; j < oridata[i].len*hashpercent; ++j)
		{
			int location = oridata[i].sig[j].location;
			for (const auto& can : hash[oridata[i].ebdqgrams[location]])
			{
				if (oridata[i].fwrv == 0 && oridata[can.id].index < oridata[i].index)
				if (edit_dp(&(oridata[i].cstr[location]), q, &(oridata[can.id].cstr[can.location]), q, 2) != -1)
					{
						if (candidate.find(can.id) == candidate.end())
						{
							vector<match> newlist;
							newlist.push_back({ location, can.location }); //(segment1,string1,segment2,string2)
							candidate[can.id] = newlist;
						}
						else
						{
							candidate[can.id].push_back({ location, can.location });
						}
					}
			}
		}

		double pretime = get_cpu_time();
		for (auto& can : candidate)
		{
			if (!DEBUG && can.second.size() >= threshold)
			{
				verify(i, can.first, can.second, DEBUG);
			}
			else if (oridata[i].index == id1 && oridata[can.first].index == id2 && oridata[can.first].fwrv == fwd)
				verify(i, can.first, can.second, DEBUG);
		}
		vtime += get_cpu_time() - pretime;
	}
}

int main(int argc, char **argv)
{
	//srand(time(NULL));
	srand(110121);

	filename = argv[1];

	double a = 0, b = 0, c = 0, d = 0;
	a = get_cpu_time();
	readfasta();
	b = get_cpu_time();
	qgramandhash();
	unordered_map<string, list<signature>> hash_table = buildtable();
	c = get_cpu_time();
	join(hash_table);
	d = get_cpu_time();

	fprintf(stderr, "Total time:%fs Initialzation time: %fs Hashing time: %fs Join time: %fs Verfication time %fs \n",d - a, b-a, c-b, d-c, vtime);

	return 0;
}