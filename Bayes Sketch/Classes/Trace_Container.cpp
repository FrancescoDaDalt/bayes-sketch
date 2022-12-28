//
//  Trace_Container.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "Trace_Container.hpp"

/*
 First constructor, which takes in a trace in form of csv file.
 */
Trace_Container::Trace_Container(string path) {
	tracetype = "Real";
	trace = path;
	set<int> id_set;
	fstream fin;
	fin.open(path, ios::in);
	string line, word, temp;
	l1 = 0;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
		id_set.insert(stoi(word));
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
		
		// Timeless data is basically the trace, but without timestamps
		timeless_data.push_back(tmp);
		
		// The Trace class computes a ground truth when reading in the file
		ground_truth[get<0>(tmp)] += get<1>(tmp);
		l1 += get<1>(tmp);
	}
	
	// The trace class also keeps a list of all distinct ids that appeared.
	for (int id : id_set) {
		distinct_ids.push_back(id);
	}
}

Trace_Container::Trace_Container(string path, int num_samples, int seed) {
	tracetype = "Real Subsampled";
	trace = path + to_string(num_samples) + " Samples";
	
	
	
	
	unordered_set<int> id_set_temp;
	vector<int> distinct_ids_temp;
	fstream fin;
	fin.open(path, ios::in);
	string line, word, temp;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
		id_set_temp.insert(stoi(word));
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
	}
	fin.close();
	
	for (int id : id_set_temp) {
		distinct_ids_temp.push_back(id);
	}
	
	shuffle(distinct_ids_temp.begin(), distinct_ids_temp.end(), default_random_engine(seed));
	int num = 0;
	unordered_set<int> id_set;
	for (int id: distinct_ids_temp) {
		id_set.insert(id);
		num ++;
		if (num == num_samples) {
			break;
		}
	}
	
	
	fin.open(path, ios::in);
	l1 = 0;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
		int id = stoi(word);
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
		if (id_set.count(id) == 1) {
			// Timeless data is basically the trace, but without timestamps
			timeless_data.push_back(tmp);
			
			// The Trace class computes a ground truth when reading in the file
			ground_truth[get<0>(tmp)] += get<1>(tmp);
			l1 += get<1>(tmp);
		}
	}
	
	// The trace class also keeps a list of all distinct ids that appeared.
	for (int id : id_set) {
		distinct_ids.push_back(id);
	}
}

Trace_Container::Trace_Container(string path, tuple<float, float> range, int bloom_filter_cells) {
	tracetype = "Real";
	trace = path;
	set<int> id_set;
	fstream fin;
	fin.open(path, ios::in);
	string line, word, temp;
	l1 = 0;
	vector<tuple<int, int>> tmp_timeless_data;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
//		id_set.insert(stoi(word));
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
		
		// Timeless data is basically the trace, but without timestamps
		tmp_timeless_data.push_back(tmp);
		
		// The Trace class computes a ground truth when reading in the file
//		ground_truth[get<0>(tmp)] += get<1>(tmp);
//		l1 += get<1>(tmp);
	}
	
	long begin = (long) (get<0>(range) * tmp_timeless_data.size());
	long end = (long) (get<1>(range) * tmp_timeless_data.size());
	timeless_data = vector<tuple<int, int>>(tmp_timeless_data.begin() + begin, tmp_timeless_data.begin() + end);
	
	for (auto tup : timeless_data) {
		id_set.insert(get<0>(tup));
		ground_truth[get<0>(tup)] += get<1>(tup);
		l1 += get<1>(tup);
	}
	
	// The trace class also keeps a list of all distinct ids that appeared.
	for (int id : id_set) {
		distinct_ids.push_back(id);
	}
	
	// Bloom filter simulation
	if (bloom_filter_cells >= 1) {
		simulateBloomFilterDistinctIDs(bloom_filter_cells);
	}
}

/*
 Second constructor which is used for syntehtic distributions.
 The arguments specifies how many items we want in the stream.
 */
Trace_Container::Trace_Container(int num_ids) {
	for (int i = 0; i < num_ids; i++) {
		distinct_ids.push_back(i);
	}
}

void Trace_Container::simulateBloomFilterDistinctIDs(int bloom_filter_cells) {
	int num_hash_functions = (int) max(1.0, ((0.7 * bloom_filter_cells) / (distinct_ids.size())));
	char* bloom_filter = (char*) calloc(sizeof(char), bloom_filter_cells);
	bloom_distinct_ids.clear();
	for (int id : distinct_ids) {
		for (int i = 0; i < num_hash_functions; i++) {
			int h = 0;
			MurmurHash3_x86_32((const void *) &id, (int) sizeof(int), (uint32_t) (i + 4224), (void *) &h);
			h = ((h % bloom_filter_cells) + bloom_filter_cells) % bloom_filter_cells;
			if (bloom_filter[h] == 0) {
				bloom_distinct_ids.push_back(id);
			}
			bloom_filter[h] = 1;
		}
	}
	free((void*) bloom_filter);
}

/*
 Synthesizes a Poisson distributed trace.
 Items are i.i.d.
 */
void Trace_Container::synthesizePoisson(int mean, int seed) {
	tracetype = "Synthetic";
	trace = "Poisson " + to_string(mean);
	timeless_data.clear();
	ground_truth.clear();
	mt19937 gen(seed);
	poisson_distribution<> d(mean);
	l1 = 0;
	for (int id : distinct_ids) {
		int sample = d(gen) + 1;
		timeless_data.push_back(make_tuple(id, sample));
		ground_truth[id] = sample;
		l1 += sample;
	}
}

void Trace_Container::synthesizeTranslatedPoisson(int mean, int translation, int seed) {
	tracetype = "Synthetic";
	trace = "Translated Poisson " + to_string(mean) + " " + to_string(translation);
	timeless_data.clear();
	ground_truth.clear();
	mt19937 gen(seed);
	poisson_distribution<> d(mean);
	l1 = 0;
	for (int id : distinct_ids) {
		int sample = (d(gen) + 1) + translation;
		timeless_data.push_back(make_tuple(id, sample));
		ground_truth[id] = sample;
		l1 += sample;
	}
}

void Trace_Container::synthesizeTranslatedGeometric(int mean, int translation, int seed) {
	tracetype = "Synthetic";
	trace = "Translated Geometric " + to_string(mean) + " " + to_string(translation);
	timeless_data.clear();
	ground_truth.clear();
	mt19937 gen(seed);
	geometric_distribution<> d(1.0 / mean);
	l1 = 0;
	for (int id : distinct_ids) {
		int sample = (d(gen) + 1) + translation;
		timeless_data.push_back(make_tuple(id, sample));
		ground_truth[id] = sample;
		l1 += sample;
	}
}

void Trace_Container::synthesizeTranslatedStudentT(int dofs, int translation, double scale, int seed) {
	tracetype = "Synthetic";
	trace = "Translated Student T " + to_string(dofs) + " " + to_string(translation) + " " + to_string(scale);
	timeless_data.clear();
	ground_truth.clear();
	mt19937 gen(seed);
	student_t_distribution<> d(dofs);
	l1 = 0;
	for (int id : distinct_ids) {
		int sample = (int) ((d(gen) + translation) * scale);
		timeless_data.push_back(make_tuple(id, sample));
		ground_truth[id] = sample;
		l1 += sample;
	}
}

void Trace_Container::synthesizeSingleOutlier(int base, int extreme) {
	tracetype = "Synthetic";
	trace = "SingleOutlier " + to_string(base) + " " + to_string(extreme);
	timeless_data.clear();
	ground_truth.clear();
	l1 = 0;
	for (int id : distinct_ids) {
		if (id == 0) {
			timeless_data.push_back(make_tuple(id, extreme));
			ground_truth[id] = extreme;
			l1 += extreme;
		} else {
			timeless_data.push_back(make_tuple(id, base));
			ground_truth[id] = base;
			l1 += base;
		}
		
	}
}

vector<tuple<int, int>> Trace_Container::get_timeless_data() {return timeless_data;}

vector<int> Trace_Container::get_distinct_ids() {return distinct_ids;}

vector<int> Trace_Container::get_filtered_distinct_ids() {
	assert(bloom_distinct_ids.size() > 0);
	return bloom_distinct_ids;
}

vector<int> Trace_Container::get_noisy_distinct_ids(double fraction) {
	auto rng = default_random_engine {};
	shuffle(begin(distinct_ids), end(distinct_ids), rng);
	vector<int> t(distinct_ids.begin(), distinct_ids.begin() + (int) (distinct_ids.size() * fraction));
	return t;
}

map<int, int> Trace_Container::get_ground_truth() {return ground_truth;}
