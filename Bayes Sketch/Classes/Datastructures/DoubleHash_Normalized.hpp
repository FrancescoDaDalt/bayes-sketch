//
//  DoubleHash_Normalized.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 26.03.22.
//

#ifndef DoubleHash_Normalized_hpp
#define DoubleHash_Normalized_hpp

#include <stdio.h>
#include <thread>
#include <armadillo>
#include "MurmurHash3.h"
#include "GlobalDatastructure_Interface.hpp"

using GDI = GlobalDatastructure_Interface;
using namespace std;

/*
 Implements the Global Datastructure interface and corresponds to the datastructure used by SeqSketch.
 */
class DoubleHash_Normalized: protected GDI  {
private:
	double* hash_tables;
	
	void inc_hash_tables(int volume, int cardinality, int hash, double factor, int level);
	void inc_vol(int volume, int hash, double factor, int level);
	void inc_card(int cardinality, int hash, int level);
	tuple<double, int> read_hash_table(int hash, int level);
	int salted_hash(int key, int level);
	
	void insert_payload(vector<tuple<int, int>> &trace, int bottom_level, int top_level, int cardinality_inc);
	void insert_payload_parallel(vector<tuple<int, int>> &trace, int cardinality_inc, int num_threads);
	void insert_ids(vector<int> &ids, int bottom_level, int top_level, int cardinality_inc);
	void insert_ids_parallel(vector<int> &ids, int cardinality_inc, int num_threads);
	double get_factor(int id, int level);
	double InverseStandardNormal(double x);
public:
	
	/*
	 Constructs a datastructure of given dimensions which normalizes entries as was proposed for the Seq-Sketch.
	 */
	DoubleHash_Normalized(int d, int w);
	
	void print();
	
	/*
	 Given a list of keys, this function returns the "sensing matrix" together with a mapping of keys to columns in the sensing matrix.
	 */
	tuple<arma::sp_dmat, map<int, int>> get_hash_matrix(vector<int> &distinct_ids);
	
	/*
	 Returns the vector b of the system Ax=b where A is the sensing matrix and x are the aggregate volumes of the keys.
	 */
	arma::dcolvec get_volume_counters();
	
	/*
	 Feeds the trace to the datastructure
	 */
	void feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise);
	
	
};

#endif /* DoubleHash_Normalized_hpp */
