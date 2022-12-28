//
//  PR_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 23.03.22.
//

#ifndef PR_Sketch_hpp
#define PR_Sketch_hpp

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include "Eigen/SparseQR"
#include "nlopt.hpp"
#include "GlobalSketch_Interface.hpp"

using GSI = GlobalSketch_Interface;
using namespace std;

class PR_Sketch : protected GSI {
private:
	
protected:
	
	static double f(unsigned int n, const double* x, double* grad, void *data);
	
	static double eq(unsigned int n, const double* x, double* grad, void *data);
	
public:
	PR_Sketch() {};
	
	/*
	 Implements the PR-Sketch method of aggregate volume estimation
	 */
	vector<double> predict(arma::sp_dmat &hash_matrix, arma::dcolvec &volume_counters);
};

#endif /* PR_Sketch_hpp */
