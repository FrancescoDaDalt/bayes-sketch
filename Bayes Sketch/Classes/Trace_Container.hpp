//
//  Trace_Container.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef Trace_Container_hpp
#define Trace_Container_hpp

#include <stdio.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <unordered_set>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include "MurmurHash3.h"

using namespace std;

/*
 A simple container class for a trace.
 Computes the list of distinct ids and the ground truth of the trace.
 */
class Trace_Container {
private:
	int l1;
	vector<tuple<int, int>> timeless_data;
	vector<int> distinct_ids;
	vector<int> bloom_distinct_ids;
	map<int, int> ground_truth;
	string tracetype;
	string trace;
	
public:
	/*
	 Read in whole trace from path
	 */
	Trace_Container(string path);
	
	/*
	 Read in a trace and randomly subsample keys from it
	 */
	Trace_Container(string path, int num_samples, int seed);
	
	/*
	 Read in a trace and take only stream elements whose rank in the stream is between x*range[0] (inclusive) and x*range[1] (exclusive) where x is the length of the trace.
	 */
	Trace_Container(string path, tuple<float, float> range, int bloom_filter_cells);
	
	/*
	 Define a trace with a number of keys. Must be followed by a synthetization step
	 */
	Trace_Container(int num_ids);
	/*
	 Synthesizes a Poisson distributed trace.
	 Items are i.i.d.
	 */
	void synthesizePoisson(int mean, int seed);
	
	/*
	 Synthesizes a translated Poisson distributed trace.
	 Items are i.i.d.
	 */
	void synthesizeTranslatedPoisson(int mean, int translation, int seed);
	
	/*
	 Synthesizes a trace where all items have size base, but one item has size extreme
	 */
	void synthesizeSingleOutlier(int base, int extreme);
	
	/*
	 Synthesizes a Tranlated Geometric trace.
	 Items are i.i.d.
	 */
	void synthesizeTranslatedGeometric(int mean, int translation, int seed);
	
	/*
	 Synthesizes a Tranlated Student-t trace.
	 Items are i.i.d.
	 */
	void synthesizeTranslatedStudentT(int dofs, int translation, double scale, int seed);
	
	/*
	 Simulates which IDs would be detected by a bloom filter of a certain size.
	 Only has false negative errors, i.e. it reports an ID as missing although it exists.
	 The number of hash functions to use is choosen to be approximately optimal,
	 given information about the number of items in the stream and the size of the filter.
	 */
	void simulateBloomFilterDistinctIDs(int bloom_filter_cells);
	
	vector<tuple<int, int>> get_timeless_data();
	vector<int> get_distinct_ids();
	vector<int> get_filtered_distinct_ids();
	vector<int> get_noisy_distinct_ids(double fraction);
	map<int, int> get_ground_truth();
	
	string get_tracetype() {return tracetype;};
	string get_trace() {return trace;};
};

#endif /* Trace_Container_hpp */
