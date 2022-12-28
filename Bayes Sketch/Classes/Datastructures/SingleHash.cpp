//
//  SingleHash.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "SingleHash.hpp"

/*
 Constructor. Requires only d and w in order to allocate the hash tables.
 */
SingleHash::SingleHash(int d, int w) {
	LDI::d = d;
	LDI::w = w;
	GDI::d = d;
	GDI::w = w;
	// The has tables are organized in such a way that volume and cardinality counter are stored next to eac other in memory.
	hash_tables = (int*) calloc(LDI::d * LDI::w * 2, sizeof(int));
}

/*
 Given an item id and a range of hash tables spcified by bottom_level and top_level, this function returns the list of
 cardinality and volume counter values associated to the given item id.
*/
vector<tuple<int, int>> SingleHash::get_vol_card(int id, int bottom_level, int top_level) {
	vector<tuple<int, int>> ret;
	for (int i = bottom_level; i < top_level; i++) {
		int hashed = salted_hash(id, i);
		ret.push_back(read_hash_table(hashed, i));
	}
	return ret;
}

/*
 Increments both cardinality and volume of a specific entry defined by the level of the hash table and the hash inside the hash table.
 Amount of increment is given by volume_inc and cardinality_inc respecitvely. If cardinality is computed by an offline-run,
 cardinality_inc should be zero in the first run when payload size is being fed to the datastrcuture.
 */
void SingleHash::inc_hash_tables(int volume_inc, int cardinality_inc, int hash, int level) {
	int index_base = 2 * (level * LDI::w + hash);
	hash_tables[index_base] += volume_inc;
	hash_tables[index_base + 1] += cardinality_inc;
}

// Same as inc_matrix but only volume increase
void SingleHash::inc_vol(int volume, int hash, int level) {
	int index_base = 2 * (level * LDI::w + hash);
	hash_tables[index_base] += volume;
}

// Same as inc_matrix but only cardinality increase
void SingleHash::inc_card(int cardinality, int hash, int level) {
	int index_base = 2 * (level * LDI::w + hash);
	hash_tables[index_base + 1] += cardinality;
}

/*
 Given the level (or depth) of the hash table and the hash inside the hash table, this function returns a tuple which contains as first entry the
 volume at the counter and second entry is the cardinality.
 */
tuple<int, int> SingleHash::read_hash_table(int hash, int level) {
	int index_base = 2 * (level * LDI::w + hash);
	return make_tuple(hash_tables[index_base], hash_tables[index_base + 1]);
}

/*
 Evaluates the hash function at a specific level, at a given key. The hash function level is used as seed to MurmurHash in order to get
 independent-ish hash functions for each level.
 */
int SingleHash::salted_hash(int key, int level) {
	int dest = 0;
	MurmurHash3_x86_32((const void *) &key, (int) sizeof(int), (uint32_t) level, (void *) &dest);
	return ((dest % LDI::w) + LDI::w) % LDI::w;
}

/*
 Takes in a trace in form of a list of tuples with the first tuple element being the flow id and the second being the payload.
 Bottom level and top level are used to define which contiguous region of hash tables should be updated. Usefull for multithreading.
 */
void SingleHash::insert_payload(vector<tuple<int, int>> &trace, int bottom_level, int top_level, int cardinality_inc) {
	for (tuple<int, int> tup : trace) {
		if (bottom_level == 0 && bottom_level != top_level) {
			LDI::l1 += get<1>(tup);
			GDI::l1 += get<1>(tup);
			l += cardinality_inc;
		}
		for (int i = bottom_level; i < top_level; i++) {
			int hashed = salted_hash(get<0>(tup), i);
			inc_hash_tables(get<1>(tup), cardinality_inc, hashed, i);
		}
	}
}

/*
 Same as insert_payload except that only cardinality counter is updated.
 Used by off-line run which computes the ground truth cardinality counters.
 */
void SingleHash::insert_ids(vector<int> &ids, int bottom_level, int top_level, int cardinality_inc) {
	for (int id : ids) {
		if (bottom_level == 0 && bottom_level != top_level) {
			l += cardinality_inc;
		}
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
void SingleHash::insert_payload_parallel(vector<tuple<int, int>> &trace, int cardinality_inc, int num_threads) {
	thread* thread_array = (thread*) calloc(num_threads, sizeof(thread));
	
	auto threadify = [this, &trace, cardinality_inc] (int bottom_level, int top_level) {insert_payload(trace, bottom_level, top_level, cardinality_inc);};
	
	for (int i = 0; i < num_threads; i++) {
		int bottom_level = i * LDI::d / num_threads;
		int top_level = (i + 1) * LDI::d / num_threads;
		thread_array[i] = thread(threadify, bottom_level, top_level);
	}
	
	for (int i = 0; i < num_threads; i++) {
		thread_array[i].join();
	}
}

/*
 Same as insert_payload_parallel except that only cardinality counter is updated.
 Used by off-line run which computes the ground truth cardinality counters.
 */
void SingleHash::insert_ids_parallel(vector<int> &ids, int cardinality_inc, int num_threads) {
	thread* thread_array = (thread*) calloc(num_threads, sizeof(thread));
	
	auto threadify = [this, &ids, cardinality_inc] (int bottom_level, int top_level) {insert_ids(ids, bottom_level, top_level, cardinality_inc);};
	
	for (int i = 0; i < num_threads; i++) {
		int bottom_level = i * LDI::d / num_threads;
		int top_level = (i + 1) * LDI::d / num_threads;
		thread_array[i] = thread(threadify, bottom_level, top_level);
	}
	
	for (int i = 0; i < num_threads; i++) {
		thread_array[i].join();
	}
}

// Takes in a trace and fully updates the whole datastructure.
void SingleHash::feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise) {
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
	insert_ids_parallel(distinct_ids, 1, 8);
}

// Prints the hash tables.
void SingleHash::print() {
	for (int i = 0; i < LDI::d; i++) {
		for (int j = 0; j < LDI::w; j++) {
			cout<<" [ "<<hash_tables[2 * (i * LDI::w + j)]<<" : "<<hash_tables[2 * (i * LDI::w + j) + 1]<<" ] \t";
		}
		cout<<endl;
	}
	cout<<endl<<LDI::l1<<endl;
}


arma::dcolvec SingleHash::get_volume_counters() {
	arma::dcolvec ret(GDI::w * GDI::d);
	for (int i = 0; i < GDI::d; i++) {
		for (int j = 0; j < GDI::w; j++) {
			ret(i * GDI::w + j) = get<0>(read_hash_table(j, i));
		}
	}
	return ret;
}


tuple<arma::sp_dmat, map<int, int>> SingleHash::get_hash_matrix(vector<int> &distinct_ids) {
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
			vals(it) = 1.0;
			it ++;
		}
		id_index_map.emplace(id, index);
		index++;
	}
	arma::sp_dmat arma_hash_matrix(locs, vals);
	return make_tuple(arma_hash_matrix, id_index_map);
}
