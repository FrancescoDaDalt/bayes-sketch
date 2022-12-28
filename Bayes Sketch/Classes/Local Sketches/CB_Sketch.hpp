//
//  CB_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef CB_Sketch_hpp
#define CB_Sketch_hpp

#include <stdio.h>
#include <limits>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "LocalSketch_Interface.hpp"
#include "LocalDatastructure_Interface.hpp"
#include "LocalSketch_Engine.hpp"

using LSI = LocalSketch_Interface;
using namespace std;

class CB_Sketch : protected LSI {
private:
	
protected:
	// Priors
	double mu;
	double chi;
	
public:
	CB_Sketch(): chi(numeric_limits<double>::infinity()), mu(0.0) {};
	
	/*
	 Grid search over all mu and chi values to minimize / maximize some statistic which is returned by the sketch engine.
	 */
	tuple<double, double> maximize_on_trace(LDI* ds, vector<int> ids, map<int, int> ground_truth, vector<double> chivalues, vector<double> muvalues, string to_optimize);
	tuple<double, double> minimize_on_trace(LDI* ds, vector<int> &ids, map<int, int> &ground_truth, vector<double> &chivalues, vector<double> &muvalues, string to_optimize);
	
	/*
	 Implements the CB method of aggregate volume estimation
	 */
	double predict(vector<tuple<int, int>> &vol_card_list);
	
	string get_prior();
	string get_priormu();
	string get_priorchi();
	
	/*
	 Sets the prior mu parameter. Default is zero
	 */
	void set_mu (double m);
	
	/*
	 Sets the prior chi parameter. Default is infinity
	 */
	void set_chi(double c);
};

#endif /* CB_Sketch_hpp */
