//
//  LocalSketch_Interface.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 16.03.22.
//

#ifndef LocalSketch_Interface_hpp
#define LocalSketch_Interface_hpp

#include <stdio.h>
#include <vector>
#include <tuple>

using namespace std;

/*
 Interface for what a local sketch algorithm (i.e. C, CB, CM, etc but not PR and the compressive sensing sketch) needs to provide.
 A local sketch needs to predict a value of the item size, given the list of volume / cardinality counters associated to that item.
 */
class LocalSketch_Interface {
private:
	
protected:
	int l;
	int l1;
	int d;
	int w;
	
public:
	LocalSketch_Interface() {};
	
	/*
	 Implementations must compute aggregate volume estimates based on a list of volume and caridnality counter measurements.
	 */
	virtual double predict(vector<tuple<int, int>> &vol_card_list) {assert(false && "Missing predict implemntation");}
	
	void set_l(int l_) {l = l_;}
	void set_l1(int l1_) {l1 = l1_;}
	void set_d(int d_) {d = d_;}
	void set_w(int w_) {w = w_;}
	
	int get_l() {return l;};
	int get_l1() {return l1;};
	int get_d() {return d;};
	int get_w() {return w;};
	
};

#endif /* LocalSketch_Interface_hpp */


