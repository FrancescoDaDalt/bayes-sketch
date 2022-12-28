//
//  C_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef C_Sketch_hpp
#define C_Sketch_hpp

#include <stdio.h>
#include <limits>
#include <cmath>
#include "LocalSketch_Interface.hpp"

using LSI = LocalSketch_Interface;
using namespace std;

class C_Sketch : protected LSI {
private:
	
protected:
	
public:
	C_Sketch() {};
	
	/*
	 Implements the Count-Sketch method of aggregate volume estimation
	 */
	double predict(vector<tuple<int, int>> &vol_card_list);
};

#endif /* C_Sketch_hpp */
