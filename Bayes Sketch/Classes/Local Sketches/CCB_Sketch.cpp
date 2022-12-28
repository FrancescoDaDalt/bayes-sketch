//
//  CCB_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "CCB_Sketch.hpp"

double CCB_Sketch::predict(vector<tuple<int, int>> &vol_card_list) {
	double d = (double) LSI::get_d();
	double l = (double) LSI::get_l();
	double l1 = (double) LSI::get_l1();
	double sum1 = mu / chi - d * l1;
	double sum2 = 0.0;
	double sum3 = 0.0;
	for (tuple<int, int> k : vol_card_list){
		double vs = (double) get<0>(k);
		double cs = (double) get<1>(k);
		sum2 += vs / (cs - 1);
		sum3 += (l - cs) / (cs - 1);
	}
	double est = (sum1 + (l - 1) * sum2) / (1 / chi + sum3);
	return est;
}

void CCB_Sketch::set_mu(double m) {mu = m;}

void CCB_Sketch::set_chi(double c) {chi = c;}

tuple<double, double> CCB_Sketch::maximize_on_trace(LDI* ds, vector<int> &ids, map<int, int> &ground_truth, vector<double> &chivalues, vector<double> &muvalues, string to_optimize) {
	
	double best = -numeric_limits<double>::infinity();
	tuple<double, double> best_mu_chi;
	
	for (double chi_prop: chivalues) {
		for (double mu_prop: muvalues) {
			cout<<"Mu: "<<mu_prop<<" Chi:"<<chi_prop<<endl;
			set_mu(mu_prop);
			set_chi(chi_prop);
			LocalSketch_Engine engine(this, ds);
			vector<tuple<string, double>>  results = engine.predict_and_benchmark(ids, ground_truth, 8, 500000);
			double score = -numeric_limits<double>::infinity();
			for (tuple<string, double> tup: results) {
				if (get<0>(tup).compare(to_optimize) == 0) {
					score = get<1>(tup);
					break;
				}
			}
			cout<<"Score: "<<score<<endl<<endl;
			if (best < score) {
				best = score;
				best_mu_chi = make_tuple(mu_prop, chi_prop);
			}
		}
	}
	set_mu(get<0>(best_mu_chi));
	set_chi(get<1>(best_mu_chi));
	return best_mu_chi;
}

tuple<double, double> CCB_Sketch::minimize_on_trace(LDI* ds, vector<int> &ids, map<int, int> &ground_truth, vector<double> &chivalues, vector<double> &muvalues, string to_optimize) {
	
	double best = numeric_limits<double>::infinity();
	tuple<double, double> best_mu_chi;
	
	for (double chi_prop: chivalues) {
		for (double mu_prop: muvalues) {
			cout<<"Mu: "<<mu_prop<<" Chi:"<<chi_prop<<endl;
			set_mu(mu_prop);
			set_chi(chi_prop);
			LocalSketch_Engine engine(this, ds);
			vector<tuple<string, double>>  results = engine.predict_and_benchmark(ids, ground_truth, 8, 500000);
			double score = numeric_limits<double>::infinity();
			for (tuple<string, double> tup: results) {
				if (get<0>(tup).compare(to_optimize) == 0) {
					score = get<1>(tup);
					break;
				}
			}
			cout<<"Score: "<<score<<endl<<endl;
			if (best > score) {
				best = score;
				best_mu_chi = make_tuple(mu_prop, chi_prop);
			}
		}
	}
	set_mu(get<0>(best_mu_chi));
	set_chi(get<1>(best_mu_chi));
	return best_mu_chi;
}

string CCB_Sketch::get_prior() {
	if (isinf(chi)) {
		return "Uninformed";
	} else {
		return ("Mu " + to_string(mu) + " Chi " + to_string(chi));
	}
}

string CCB_Sketch::get_priormu() {
	stringstream stream;
	stream<<fixed<<setprecision(20)<<mu;
	string s = stream.str();
	return s;
}

string CCB_Sketch::get_priorchi() {
	stringstream stream;
	stream<<fixed<<setprecision(20)<<chi;
	string s = stream.str();
	return s;
}
