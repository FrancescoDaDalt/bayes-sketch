//
//  SingleHash.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef SingleHash_hpp
#define SingleHash_hpp

#include <stdio.h>
#include <thread>
#include <armadillo>
#include "MurmurHash3.h"
#include "LocalDatastructure_Interface.hpp"
#include "GlobalDatastructure_Interface.hpp"

using LDI = LocalDatastructure_Interface;
using GDI = GlobalDatastructure_Interface;
using namespace std;

/*
 Implements the Local and Gloval Datastructure interface and corresponds to the datastructure used by CM, CB, CCB and CCA sketch, and by PR and SeqSketch.
 */
class SingleHash: protected LDI, protected GDI {
private:
	int* hash_tables;
	
	void inc_hash_tables(int volume, int cardinality, int hash, int level);
	void inc_vol(int volume, int hash, int level);
	void inc_card(int cardinality, int hash, int level);
	tuple<int, int> read_hash_table(int hash, int level);
	int salted_hash(int key, int level);
	
	void insert_payload(vector<tuple<int, int>> &trace, int bottom_level, int top_level, int cardinality_inc);
	void insert_payload_parallel(vector<tuple<int, int>> &trace, int cardinality_inc, int num_threads);
	void insert_ids(vector<int> &ids, int bottom_level, int top_level, int cardinality_inc);
	void insert_ids_parallel(vector<int> &ids, int cardinality_inc, int num_threads);
public:
	using LDI::get_l1;
	using LDI::get_l;
	/*
	 Constructs a datastructure of given dimensions with no special behavior. This is the datastructure used by CountMin, CB, CCB, CCA and PR Sketch.
	 */
	SingleHash(int d, int w);
	void print();
	
	/*
	 Returns the vector b of the system Ax=b where A is the sensing matrix and x are the aggregate volumes of the keys.
	 */
	vector<tuple<int, int>> get_vol_card(int id, int bottom_level, int top_level);
	
	/*
	 Given a list of keys, this function returns the "sensing matrix" together with a mapping of keys to columns in the sensing matrix.
	 */
	tuple<arma::sp_dmat, map<int, int>> get_hash_matrix(vector<int> &distinct_ids);
	
	/*
	 Returns a list of <volume | cardinality> tuples associated to a certain key.
	 Tuples are collected for datastructure counter arrays between bottom and top level.
	 */
	arma::dcolvec get_volume_counters();
	
	/*
	 Feeds the trace to the datastructure
	 */
	void feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise);
	
	
};

#endif /* SingleHash_hpp */
