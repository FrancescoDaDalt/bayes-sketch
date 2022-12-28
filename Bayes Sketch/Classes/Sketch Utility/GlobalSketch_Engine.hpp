//
//  GlobalSketch_Engine.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 22.03.22.
//

#ifndef GlobalSketch_Engine_hpp
#define GlobalSketch_Engine_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <thread>
#include <numeric>
#include <future>
#include "Eigen/Dense"

#include "GlobalSketch_Interface.hpp"
#include "GlobalDatastructure_Interface.hpp"

using GSI = GlobalSketch_Interface;
using GDI = GlobalDatastructure_Interface;
using namespace std;

/*
 The idea is to combine bolier-plate code for actually executing the sketch on a datastrcuture.
 Provides a common test-bed for all sketches.
 */
class GlobalSketch_Engine {
private:
	GSI* sketch;
	GDI* ds;
	
protected:
	vector<tuple<int, vector<tuple<int, int>>>> collect_data(vector<int> ids, int bottom_level, int top_level);
	
public:
	
	/*
	 Constructs the engine which requires a global sketch object and a global datastructure object
	 */
	GlobalSketch_Engine(GSI* sketch, GDI* ds);
	
	/*
	 Predicts, based on the sketch and datastructure, the aggregate volumes of the ids passed as argument and computes error statistics based in the ground truth passed as parameter.
	 */
	vector<tuple<string, double>> predict_and_benchmark(vector<int> ids, map<int, int> ground_truth);
	
};

#endif /* GlobalSketch_Engine_hpp */
