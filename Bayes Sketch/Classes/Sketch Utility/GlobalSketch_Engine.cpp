//
//  GlobalSketch_Engine.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 22.03.22.
//

#include "GlobalSketch_Engine.hpp"

// Constructor. Handles setting the right 'static' parameters of the sketch algorithm.
GlobalSketch_Engine::GlobalSketch_Engine(GSI* sketch, GDI* ds): sketch(sketch), ds(ds) {
	sketch->set_d(ds->get_d());
	sketch->set_w(ds->get_w());
	sketch->set_l1(ds->get_l1());
}

vector<tuple<string, double>> GlobalSketch_Engine::predict_and_benchmark(vector<int> ids, map<int, int> ground_truth) {
	
	
	auto timer_start = chrono::steady_clock::now();
	int l1 = ds->get_l1();
	
	auto hash_matrix_and_id_index_map = ds->get_hash_matrix(ids);
	arma::sp_dmat hash_matrix = get<0>(hash_matrix_and_id_index_map);
	map<int, int> id_index_map = get<1>(hash_matrix_and_id_index_map);

	arma::dcolvec volume_counters = ds->get_volume_counters();

	vector<double> ans = sketch->predict(hash_matrix, volume_counters);
	
	Eigen::MatrixXd gt_est_array(ids.size(), 2);
	
	int global_index = 0;
	for (int id: ids) {
		gt_est_array(global_index, 0) = (double) ground_truth.at(id);
		gt_est_array(global_index, 1) = ans.at(id_index_map.at(id));
		global_index++;
	}
	auto timer_end = chrono::steady_clock::now();

	vector<tuple<string, double>> return_statistics;
	
	Eigen::MatrixXd abs_err_array = (gt_est_array.col(0) - gt_est_array.col(1)).array().abs().matrix();
	Eigen::MatrixXd err_array = (gt_est_array.col(0) - gt_est_array.col(1));
	double avg_abs_error = abs_err_array.array().mean();
	
	double avg_eps_error = abs_err_array.array().mean() / l1;
	
	auto median_it = abs_err_array.data() + abs_err_array.size() / 2;
	nth_element(abs_err_array.data(), median_it , abs_err_array.data() + abs_err_array.size());
	double med_eps_error = (*median_it) / l1;
	
	Eigen::MatrixXd rel_err_array = (((gt_est_array.col(0) - gt_est_array.col(1)).array()) / gt_est_array.col(0).array()).abs().matrix();
	double avg_rel_error = rel_err_array.array().mean();
	
	median_it = rel_err_array.data() + rel_err_array.size() / 2;
	nth_element(rel_err_array.data(), median_it , rel_err_array.data() + rel_err_array.size());
	double med_rel_error = (*median_it);
	
	
	double norm_rel_error = ((gt_est_array.col(0) - gt_est_array.col(1)).norm()) / gt_est_array.col(0).norm();
	
	Eigen::MatrixXd squared_err_array = ((gt_est_array.col(0) - gt_est_array.col(1)).array()).square().matrix();
	double avg_sq_error = squared_err_array.array().mean();
	
	median_it = squared_err_array.data() + squared_err_array.size() / 2;
	nth_element(squared_err_array.data(), median_it , squared_err_array.data() + squared_err_array.size());
	double med_sq_error = (*median_it);
	
	Eigen::MatrixXd centered_gt_array = (gt_est_array.col(0).array() - gt_est_array.col(0).array().mean()).matrix();
	Eigen::MatrixXd centered_est_array = (gt_est_array.col(1).array() - gt_est_array.col(1).array().mean()).matrix();
	double cross_corr = (centered_gt_array.array() * centered_est_array.array()).sum() / (centered_gt_array.norm() * centered_est_array.norm());
	
	double rsquare_surrogate = (gt_est_array.col(0) - gt_est_array.col(1)).array().square().sum() / centered_gt_array.squaredNorm();
	
	double explained_var_surrogate = (err_array.array() - err_array.array().mean()).square().mean() / centered_gt_array.array().square().mean();
	
	Eigen::MatrixXd gt_array = gt_est_array.col(0);
	for (int percentile = 1; percentile < 100; percentile *= 2) {
		auto percentile_it = gt_array.data() + (gt_array.size() * (100 - percentile)) / 100;
		nth_element(gt_array.data(), percentile_it , gt_array.data() + gt_array.size());
		double gt_percentile = (*percentile_it);
		vector<double> tmp;
		double avg_perc_rel_error = 0.0;
		int counter = 0;
		for (int row = 0; row < gt_est_array.rows(); row ++) {
			if (gt_est_array(row, 0) >= gt_percentile) {
				avg_perc_rel_error += abs((gt_est_array(row, 0) - gt_est_array(row, 1)) / gt_est_array(row, 0));
				counter ++;
			}
		}
		avg_perc_rel_error /= counter;
		
		string tmpstr1 = "avg_rel_error " + to_string(percentile) + "'th percentile";
		string tmpstr2 = to_string(percentile) + "'th percentile";
		
		//		return_statistics.push_back(make_tuple("percentile value", (double) gt_percentile));
		return_statistics.push_back(make_tuple(tmpstr1, avg_perc_rel_error));
		return_statistics.push_back(make_tuple(tmpstr2,  gt_percentile));
		
	}
	
	
	
	
	return_statistics.push_back(make_tuple("query_time", (double) chrono::duration_cast<chrono::nanoseconds>(timer_end - timer_start).count()));
	return_statistics.push_back(make_tuple("avg_abs_error", avg_abs_error));
	return_statistics.push_back(make_tuple("avg_eps_error", avg_eps_error));
	return_statistics.push_back(make_tuple("med_eps_error", med_eps_error));
	return_statistics.push_back(make_tuple("avg_rel_error", avg_rel_error));
	return_statistics.push_back(make_tuple("med_rel_error", med_rel_error));
	return_statistics.push_back(make_tuple("avg_sq_error", avg_sq_error));
	return_statistics.push_back(make_tuple("med_sq_error", med_sq_error));
	return_statistics.push_back(make_tuple("crss_corr", cross_corr));
	return_statistics.push_back(make_tuple("norm_rel_error", norm_rel_error));
	return_statistics.push_back(make_tuple("rsquare_surrogate", rsquare_surrogate));
	return_statistics.push_back(make_tuple("explained_var_surrogate", explained_var_surrogate));
	
	return return_statistics;
}
