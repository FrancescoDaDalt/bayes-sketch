//
//  LocalSketch_Engine.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef LocalSketch_Engine_hpp
#define LocalSketch_Engine_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <thread>
#include <numeric>
#include <chrono>
#include <future>
#include "Eigen/Dense"

#include "LocalSketch_Interface.hpp"
#include "LocalDatastructure_Interface.hpp"

using LSI = LocalSketch_Interface;
using LDI = LocalDatastructure_Interface;
using namespace std;

/*
 The idea is to combine bolier-plate code for actually executing the sketch on a datastrcuture.
 Provides a common test-bed for all sketches.
 Paralellization right now is not very good because the porgram only uses multithreading for evaluating all the hashes and collecting the volume / cardinaity counter values associated to each item. The prediction per se is still done sequentially.
 */
class LocalSketch_Engine {
private:
	LSI* sketch;
	LDI* ds;
	
protected:
	vector<tuple<int, vector<tuple<int, int>>>> collect_data(vector<int> ids, int bottom_level, int top_level);
	
public:
	/*
	 Constructs the engine which requires a local sketch object and a local datastructure object
	 */
	LocalSketch_Engine(LSI* sketch, LDI* ds);
	
	/*
	 Predicts, based on the sketch and datastructure, the aggregate volumes of the ids passed as argument and computes error statistics based in the ground truth passed as parameter.
	 Makes use of parallelization and batching although the implementation is suboptimal for huge traces in all honesty.
	 */
	vector<tuple<string, double>> predict_and_benchmark(vector<int> ids, map<int, int> ground_truth, int num_threads, int batch_size);
	
};

#endif /* LocalSketch_Engine_hpp */
