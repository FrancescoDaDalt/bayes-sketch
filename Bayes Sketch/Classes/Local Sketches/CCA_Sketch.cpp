//
//  CCA_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "CCA_Sketch.hpp"

double CCA_Sketch::predict(vector<tuple<int, int>> &vol_card_list) {
	double w = (double) LSI::get_w();
	double d = (double) LSI::get_d();
	double l = (double) LSI::get_l();
	double l1 = (double) LSI::get_l1();
	double sum = 0.0;
	for (tuple<int, int> k : vol_card_list){
		double vs = (double) get<0>(k);
		double cs = (double) get<1>(k);
		sum += ((vs / cs) * (l + w - 1) - l1) / (w - 1);
	}
	double est = sum / d;
	return est;
}
