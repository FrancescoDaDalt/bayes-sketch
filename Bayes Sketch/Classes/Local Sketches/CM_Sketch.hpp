//
//  CM_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef CM_Sketch_hpp
#define CM_Sketch_hpp

#include <stdio.h>
#include <limits>
#include <cmath>
#include "LocalSketch_Interface.hpp"

using LSI = LocalSketch_Interface;
using namespace std;

class CM_Sketch : protected LSI {
private:
	
protected:
	
public:
	CM_Sketch() {};
	
	/*
	 Implements the CountMin method of aggregate volume estimation
	 */
	double predict(vector<tuple<int, int>> &vol_card_list);
};

#endif /* CM_Sketch_hpp */
