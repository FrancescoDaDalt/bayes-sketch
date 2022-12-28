//
//  LocalDatastructure_Interface.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef LocalDatastructure_Interface_hpp
#define LocalDatastructure_Interface_hpp

#include <stdio.h>
#include <vector>
#include <tuple>
#include "Trace_Container.hpp"

using namespace std;

/*
 Interface that should be implemented by datastructures that work with 'local' sketches such as C, CM, CCA, etc but not PR and compressive sensing ones.
 */
class LocalDatastructure_Interface {
protected:
	int d = 0;
	int w = 0;
	int l1 = 0;
	int l = 0;
	
	
public:
	/*
	 Returns the vector b of the system Ax=b where A is the sensing matrix and x are the aggregate volumes of the keys.
	 */
	virtual vector<tuple<int, int>> get_vol_card(int id, int bottom_level, int top_level) {assert(false && "Missing implementation of vol_card retrival");};
	
	/*
	 Feeds the trace to the datastructure
	 */
	virtual void feed_trace(Trace_Container* trace, bool use_bloom_filter, float cardinality_noise) {assert(false && "Missing implementation of trace feeding");};
	
	int get_d();
	int get_w();
	int get_l1();
	int get_l();
	
};

#endif /* LocalDatastructure_Interface_hpp */
