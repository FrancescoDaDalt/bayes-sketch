//
//  CCA_Sketch.hpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#ifndef CCA_Sketch_hpp
#define CCA_Sketch_hpp

#include <stdio.h>
#include <limits>
#include "LocalSketch_Interface.hpp"

using LSI = LocalSketch_Interface;
using namespace std;

class CCA_Sketch : protected LSI {
private:
	
protected:
	
public:
	CCA_Sketch() {};
	/*
	 Implements the CCA method of aggregate volume estimation
	 */
	double predict(vector<tuple<int, int>> &vol_card_list);
};

#endif /* CCA_Sketch_hpp */
