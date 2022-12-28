//
//  Seq_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 24.03.22.
//

#ifndef Seq_Sketch_hpp
#define Seq_Sketch_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "Eigen/SparseQR"
#include "nlopt.hpp"
#include "GlobalSketch_Interface.hpp"
#include <armadillo>
#include <KL1pInclude.h>

using GSI = GlobalSketch_Interface;
using namespace std;

class Seq_Sketch : protected GSI {
private:
	
protected:
	
	static double f(unsigned int n, const double* x, double* grad, void *data);
	
	static double eq(unsigned int n, const double* x, double* grad, void *data);
	
public:
	Seq_Sketch() {};
	
	/*
	 Implements the Seq-Sketch method of aggregate volume estimation
	 */
	vector<double> predict(arma::sp_dmat &hash_matrix, arma::dcolvec &volume_counters);
};

#endif /* Seq_Sketch_hpp */
