//
//  C_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "C_Sketch.hpp"
#include <functional>

double C_Sketch::predict(vector<tuple<int, int>> &vol_card_list) {
	vector<int> t;
	for (tuple<int, int> k : vol_card_list){
		t.push_back(get<0>(k));
	}
	int median_index = LSI::get_d() / 2;
	nth_element(t.begin(), t.begin() + median_index, t.end());
	double est = (double) t[median_index];
	return est;
}


