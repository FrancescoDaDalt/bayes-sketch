//
//  DoubleHash_Normalized.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 26.03.22.
//

#include "DoubleHash_Normalized.hpp"

/*
 Constructor. Requires only d and w in order to allocate the hash tables.
 */
DoubleHash_Normalized::DoubleHash_Normalized(int d, int w) {
	GDI::d = d;
	GDI::w = w;
	// The has tables are organized in such a way that volume and cardinality counter are stored next to eac other in memory.
	hash_tables = (double*) calloc(GDI::d * GDI::w * 2, sizeof(double));
}

/*
 Increments both cardinality and volume of a specific entry defined by the level of the hash table and the hash inside the hash table.
 Amount of increment is given by volume_inc and cardinality_inc respecitvely. If cardinality is computed by an offline-run,
 cardinality_inc should be zero in the first run when payload size is being fed to the datastrcuture.
 factor handles the second level hash.
 */
void DoubleHash_Normalized::inc_hash_tables(int volume_inc, int cardinality_inc, int hash, double factor, int level) {
	int index_base = 2 * (level * GDI::w + hash);
	hash_tables[index_base] += volume_inc * factor;
	hash_tables[index_base + 1] += cardinality_inc;
}

// Same as inc_matrix but only volume increase
void DoubleHash_Normalized::inc_vol(int volume, int hash, double factor, int level) {
	int index_base = 2 * (level * GDI::w + hash);
	hash_tables[index_base] += volume * factor;
}

// Same as inc_matrix but only cardinality increase
void DoubleHash_Normalized::inc_card(int cardinality, int hash, int level) {
	int index_base = 2 * (level * GDI::w + hash);
	hash_tables[index_base + 1] += cardinality;
}

/*
 Given the level (or depth) of the hash table and the hash inside the hash table, this function returns a tuple which contains as first entry the
 volume at the counter and second entry is the cardinality.
 */
tuple<double, int> DoubleHash_Normalized::read_hash_table(int hash, int level) {
	int index_base = 2 * (level * GDI::w + hash);
	return make_tuple(hash_tables[index_base], hash_tables[index_base + 1]);
}

/*
 Evaluates the hash function at a specific level, at a given key. The hash function level is used as seed to MurmurHash in order to get
 independent-ish hash functions for each level.
 Also used for second level hashing.
 */
int DoubleHash_Normalized::salted_hash(int key, int level) {
	int dest = 0;
	MurmurHash3_x86_32((const void *) &key, (int) sizeof(int), (uint32_t) level, (void *) &dest);
	return ((dest % GDI::w) + GDI::w) % GDI::w;
}


// From here: https://www.johndcook.com/blog/cpp_phi_inverse/
double DoubleHash_Normalized::InverseStandardNormal(double x) {
	double c[] = {2.515517, 0.802853, 0.010328};
	double d[] = {1.432788, 0.189269, 0.001308};
	
	auto RationalApproximation = [c, d] (double x) {return x - ((c[2] * x + c[1]) * x + c[0]) / (((d[2] * x + d[1]) * x + d[0]) * x + 1.0);};
	
	if (x < 0.5) {
		return -RationalApproximation( sqrt(-2.0 * log(x)));
	} else {
		return RationalApproximation( sqrt(-2.0 * log(1-x)));
	}
}

double DoubleHash_Normalized::get_factor(int id, int level) {
	int dest = 0;
	MurmurHash3_x86_32((const void *) &id, (int) sizeof(int), (uint32_t) level, (void *) &dest);
	double samplepoint = (((dest % 10000) + 10000) % 10000 + 0.5) / 10000.0;
	double standardnormalsample = InverseStandardNormal(samplepoint);
	double transformed_sample = standardnormalsample / (5.0 * sqrt(GDI::d)) + 1.0 / sqrt(GDI::d);
	return transformed_sample;
}

/*
 Takes in a trace in form of a list of tuples with the first tuple element being the flow id and the second being the payload.
 Bottom level and top level are used to define which contiguous region of hash tables should be updated. Usefull for multithreading.
 factor accounts for the second level hashing.
 */
void DoubleHash_Normalized::insert_payload(vector<tuple<int, int>> &trace, int bottom_level, int top_level, int cardinality_inc) {
	for (tuple<int, int> tup : trace) {
		if (bottom_level == 0 && bottom_level != top_level) {
			l1 += get<1>(tup);
		}
		for (int i = bottom_level; i < top_level; i++) {
			int hashed = salted_hash(get<0>(tup), i);
			double factor = get_factor(get<0>(tup), i + GDI::d);
			inc_hash_tables(get<1>(tup), cardinality_inc, hashed, factor, i);
		}
	}
}

/*
 Same as insert_payload except that only cardinality counter is updated.
 Used by off-line run which computes the ground truth cardinality counters.
 Not really needed for the Count Sketch.
 */
void DoubleHash_Normalized::insert_ids(vector<int> &ids, int bottom_level, int top_level, int cardinality_inc) {
	for (int id : ids) {
		for (int i = bottom_level; i < top_level; i++) {
			int hashed = salted_hash(id, i);
			inc_card(cardinality_inc, hashed, i);
		}
	}
}

/*
 Takes in a trace in form of list of tuples and feeds the trace to the datastructure parallelized over num_threads many threads.
 Parallelization speedup should be almost linear.
 */
void DoubleHash_Normalized::insert_payload_parallel(vector<tuple<int, int>> &trace, int cardinality_inc, int num_threads) {
	thread* thread_array = (thread*) calloc(num_threads, sizeof(thread));
	
	auto threadify = [this, &trace, cardinality_inc] (int bottom_level, int top_level) {insert_payload(trace, bottom_level, top_level, cardinality_inc);};
	
	for (int i = 0; i < num_threads; i++) {
		int bottom_level = i * d / num_threads;
		int top_level = (i + 1) * d / num_threads;
		thread_array[i] = thread(threadify, bottom_level, top_level);
	}
	
	for (int i = 0; i < num_threads; i++) {
		thread_array[i].join();
	}
}

/*
 Same as insert_payload_parallel except that only cardinality counter is updated.
 Used by off-line run which computes the ground truth cardinality counters.
 Not really needed for the Count Sketch.
 */
void DoubleHash_Normalized::insert_ids_parallel(vector<int> &ids, int cardinality_inc, int num_threads) {
	thread* thread_array = (thread*) calloc(num_threads, sizeof(thread));
	
	auto threadify = [this, &ids, cardinality_inc] (int bottom_level, int top_level) {insert_ids(ids, bottom_level, top_level, cardinality_inc);};
	
	for (int i = 0; i < num_threads; i++) {
		int bottom_level = i * d / num_threads;
		int top_level = (i + 1) * d / num_threads;
		thread_array[i] = thread(threadify, bottom_level, top_level);
	}
	
	for (int i = 0; i < num_threads; i++) {
		thread_array[i].join();
	}
}


// Takes in a trace and fully updates the whole datastructure.
void DoubleHash_Normalized::feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise) {
	vector<tuple<int, int>> timeless_data = trace->get_timeless_data();
	vector<int> distinct_ids;
	if (use_bloom_filter) {
		distinct_ids = trace->get_filtered_distinct_ids();
	} else {
		if (cardinality_noise == 0.0) {
			distinct_ids = trace->get_distinct_ids();
		} else {
			distinct_ids = trace->get_noisy_distinct_ids(1.0 - cardinality_noise);
		}
	}
	map<int, int> ground_truth = trace->get_ground_truth();
	
	insert_payload_parallel(timeless_data, 0, 8);
	
	// Can in theory be omitted.
	insert_ids_parallel(distinct_ids, 1, 8);
	
}

// Prints the hash tables.
void DoubleHash_Normalized::print() {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < w; j++) {
			cout<<" [ "<<hash_tables[2 * (i * w + j)]<<" : "<<hash_tables[2 * (i * w + j) + 1]<<" ] \t";
		}
		cout<<endl;
	}
	cout<<endl<<l1<<endl;
}


arma::dcolvec DoubleHash_Normalized::get_volume_counters() {
	arma::dcolvec ret(GDI::w * GDI::d);
	for (int i = 0; i < GDI::d; i++) {
		for (int j = 0; j < GDI::w; j++) {
			ret(i * GDI::w + j) = get<0>(read_hash_table(j, i));
		}
	}
	return ret;
}

tuple<arma::sp_dmat, map<int, int>> DoubleHash_Normalized::get_hash_matrix(vector<int> &distinct_ids) {
	map<int, int> id_index_map;
	
	arma::umat locs(2, GDI::d * distinct_ids.size());
	arma::dcolvec vals(GDI::d * distinct_ids.size());
	int index = 0;
	int it = 0;
	for (int id : distinct_ids) {
		for (int j = 0; j < GDI::d; j++) {
			int global_index = GDI::w * j + salted_hash(id, j);
			locs(0, it) = global_index;
			locs(1, it) = index;
			vals(it) = get_factor(id, j + GDI::d);
			it ++;
		}
		id_index_map.emplace(id, index);
		index++;
	}
	arma::sp_dmat arma_hash_matrix(locs, vals);
	return make_tuple(arma_hash_matrix, id_index_map);
}
