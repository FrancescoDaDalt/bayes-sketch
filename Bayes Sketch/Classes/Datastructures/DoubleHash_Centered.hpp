//
//  DoubleHash.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 18.03.22.
//

#ifndef DoubleHash_hpp
#define DoubleHash_hpp

#include <stdio.h>
#include <thread>
#include "MurmurHash3.h"
#include "LocalDatastructure_Interface.hpp"

using LDI = LocalDatastructure_Interface;
using namespace std;

/*
 Implements the Local Datastructure interface and corresponds to the datastructure used by the Count Sketch.
 Is very similar to SingleHash but has slightly different insertion and extraction formulae.
 */
class DoubleHash_Centered: protected LDI  {
private:
	int* hash_tables;
	
	void inc_hash_tables(int volume, int cardinality, int hash, int factor, int level);
	void inc_vol(int volume, int hash, int factor, int level);
	void inc_card(int cardinality, int hash, int level);
	tuple<int, int> read_matrix(int hash, int factor, int level);
	int salted_hash(int key, int level);
	
	void insert_payload(vector<tuple<int, int>> &trace, int bottom_level, int top_level, int cardinality_inc);
	void insert_payload_parallel(vector<tuple<int, int>> &trace, int cardinality_inc, int num_threads);
	void insert_ids(vector<int> &ids, int bottom_level, int top_level, int cardinality_inc);
	void insert_ids_parallel(vector<int> &ids, int cardinality_inc, int num_threads);
public:
	/*
	 Constructs a datastructure of given dimensions which multiplies the entries by +-1 as is done by the Count-Sketch
	 */
	DoubleHash_Centered(int d, int w);
	void print();
	
	/*
	 Returns a list of <volume | cardinality> tuples associated to a certain key.
	 Tuples are collected for datastructure counter arrays between bottom and top level.
	 Multiplication by +-1 is already taken care of.
	 */
	vector<tuple<int, int>> get_vol_card(int id, int bottom_level, int top_level);
	
	/*
	 Feeds the trace to the datastructure
	 */
	void feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise);
	
	
};

#endif /* DoubleHash_hpp */
