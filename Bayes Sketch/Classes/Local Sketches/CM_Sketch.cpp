//
//  CM_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "CM_Sketch.hpp"

double CM_Sketch::predict(vector<tuple<int, int>> &vol_card_list) {
	double est = numeric_limits<double>::infinity();
	for (tuple<int, int> k : vol_card_list){
		est = min(est, (double) get<0>(k));
	}
	return est;
}
