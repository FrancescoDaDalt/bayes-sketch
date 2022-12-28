//
//  GlobalDatastructure_Interface.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 22.03.22.
//

#ifndef GlobalDatastructure_Interface_hpp
#define GlobalDatastructure_Interface_hpp

#include <stdio.h>
#include <vector>
#include <tuple>
#include <armadillo>
#include "Trace_Container.hpp"
#include "Eigen/Sparse"

using namespace std;

/*
 Interface that should be implemented by datastructures that work with 'local' sketches such as C, CM, CCA, etc but not PR and compressive sensing ones.
 */
class GlobalDatastructure_Interface {
protected:
	int d = 0;
	int w = 0;
	int l1 = 0;
	
	
public:
	/*
	 Given a list of keys, this function returns the "sensing matrix" together with a mapping of keys to columns in the sensing matrix.
	 */
	virtual tuple<arma::sp_dmat, map<int, int>> get_hash_matrix(vector<int> &distinct_ids) {assert(false && "Missing implementation of hash_matrix retrival");};
	
	/*
	 Returns the vector b of the system Ax=b where A is the sensing matrix and x are the aggregate volumes of the keys.
	 */
	virtual arma::dcolvec get_volume_counters() {assert(false && "Missing implementation of volume_counter retrival");};
	
	/*
	 Feeds the trace to the datastructure
	 */
	virtual void feed_trace(Trace_Container* trace, bool use_bloom_filter) {assert(false && "Missing implementation of trace feeding");};
	
	int get_d();
	int get_w();
	int get_l1();
	
};

#endif /* GlobalDatastructure_Interface_hpp */
