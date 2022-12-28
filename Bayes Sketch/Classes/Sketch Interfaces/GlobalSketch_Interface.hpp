//
//  GlobalSketch_Interface.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 22.03.22.
//

#ifndef GlobalSketch_Interface_hpp
#define GlobalSketch_Interface_hpp

#include <stdio.h>
#include <armadillo>
#include "Eigen/Sparse"
#include "Eigen/Dense"

using namespace std;

/*
 Interface for what a global sketch algorithm (i.e. PR and compressive sensing sketch) needs to provide.
 */
class GlobalSketch_Interface {
private:
	
protected:
	int l1;
	int d;
	int w;
	
public:
	GlobalSketch_Interface() {};
	
	/*
	 Implementations must compute aggregate volume estimates based on the sensing matrix and its system RHS vector
	 */
	virtual vector<double> predict(arma::sp_dmat &hash_matrix, arma::dcolvec &volume_counters) {assert(false && "Missing predict implemntation");}
	
	void set_d(int d_) {d = d_;}
	void set_w(int w_) {w = w_;}
	void set_l1(int l1_) {l1 = l1_;}
	
	int get_d() {return d;};
	int get_w() {return w;};
	int get_l1() {return l1;};
	
};

#endif /* GlobalSketch_Interface_hpp */
