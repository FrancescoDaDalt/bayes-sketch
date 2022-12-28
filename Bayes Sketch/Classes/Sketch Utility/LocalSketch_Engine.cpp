//
//  LocalSketch_Engine.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "LocalSketch_Engine.hpp"

// Constructor. Handles setting the right 'static' parameters of teh sketch algorithm.
LocalSketch_Engine::LocalSketch_Engine(LSI* sketch, LDI* ds): sketch(sketch), ds(ds) {
	sketch->set_l(ds->get_l());
	sketch->set_l1(ds->get_l1());
	sketch->set_d(ds->get_d());
	sketch->set_w(ds->get_w());
}

/*
 Computes predictions and benchmarks them immediately. Batch size is used to control how many items should be predicted at a time (typically 100k to 1M).
 Parallelization is not very good because right now only the evaluation of hash functions is parallelized to extract for each item in the batch, the list of its associated volume and cardinality counters for a given range of hash fucntion. Each threads handles a different range of hash fucntions. After parallel execution, the volume and cardianlity counter values returned by each thread are stitched together for each item.
 The prediction given the list of volume and cardinality counter values then happens sequentially.
 */
vector<tuple<string, double>> LocalSketch_Engine::predict_and_benchmark(vector<int> ids, map<int, int> ground_truth, int num_threads, int batch_size) {
	
	auto timer_start = chrono::steady_clock::now();
	thread* thread_array = (thread*) calloc(num_threads, sizeof(thread));
	
	future<vector<tuple<int, vector<tuple<int, int>>>>>* thread_return_array = (future<vector<tuple<int, vector<tuple<int, int>>>>>*) calloc(num_threads, sizeof(future<vector<tuple<int, vector<tuple<int, int>>>>>));
	
	int cursor = min(batch_size, (int) ids.size());
	vector<int>::const_iterator bottom = ids.begin();
	vector<int>::const_iterator top = ids.begin() + cursor;
	
	Eigen::MatrixXd gt_est_array(ids.size(), 2);
	
	int d = ds->get_d();
	int l1 = ds->get_l1();
	
	int it = 0;
	
	while (bottom != top) {
		
		vector<int> id_batch(bottom, top);
		bottom = top;
		cursor = min(cursor + batch_size, (int) ids.size());
		top = ids.begin() + cursor;
		
		vector<tuple<int, vector<tuple<int, int>>>> aggregator;
		
		auto threadify = [this, id_batch] (promise<vector<tuple<int, vector<tuple<int, int>>>>> && prom, int bottom_level, int top_level) {prom.set_value(collect_data(id_batch, bottom_level, top_level));};
		
		for (int i = 0; i < num_threads; i++) {
			int bottom_level = i * d / num_threads;
			int top_level = (i + 1) * d / num_threads;
			auto a = promise<vector<tuple<int, vector<tuple<int, int>>>>>();
			thread_return_array[i] = a.get_future();
			thread_array[i] = thread(threadify, move(a), bottom_level, top_level);
		}
		
		vector<vector<tuple<int, vector<tuple<int, int>>>>> temp_aggregator;
		temp_aggregator.reserve(num_threads);
		
		for (int i = 0; i < num_threads; i++) {
			thread_array[i].join();
			temp_aggregator.push_back(thread_return_array[i].get());
		}
		
		// Stitches goether returned lists by threads.
		aggregator = temp_aggregator[0];
		for (int i = 0; i < temp_aggregator[0].size(); i++) {
			get<1>(aggregator[i]).reserve(d);
			for (int j = 1; j < num_threads; j++) {
				assert(get<0>(aggregator[i]) == get<0>(temp_aggregator[j][i]));
				get<1>(aggregator[i]).insert(get<1>(aggregator[i]).end(), get<1>(temp_aggregator[j][i]).begin(), get<1>(temp_aggregator[j][i]).end());
			}
		}
		
		for (tuple<int, vector<tuple<int, int>>> agg : aggregator) {
			// Sketch prediction happens here
			double est = sketch->predict(get<1>(agg));
			gt_est_array(it, 0) = (double) ground_truth.at(get<0>(agg));
			gt_est_array(it, 1) = est;
			it++;
//			gt_pred_list.push_back(make_tuple((double) ground_truth.at(get<0>(agg)), est));
//			error_list.push_back(ground_truth.at(get<0>(agg)) - est);
//			eps_error += abs(est - ground_truth.at(get<0>(agg))) / l1;
//			rel_error += abs(est - ground_truth.at(get<0>(agg))) / ground_truth.at(get<0>(agg));
		}
	}
	
	free(thread_array);
	free(thread_return_array);
	
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
	
	
//	cout<<"l2 = "<<gt_est_array.col(0).array().square().sum()<<endl;
//	cout<<"l1 = "<<gt_est_array.col(0).array().sum()<<endl;
//	cout<<"l1 = "<<l1<<endl;
	
	return return_statistics;
}

/*
 Returns for a list of ids, a list which contains lists which indicate the volume and cardinality counters associated to the ids, where the counters are evaluated only in a specified range of hash functions.
 */
vector<tuple<int, vector<tuple<int, int>>>> LocalSketch_Engine::collect_data(vector<int> ids, int bottom_level, int top_level) {
	vector<tuple<int, vector<tuple<int, int>>>> ret;
	for (int id : ids) {
		ret.insert(ret.end(), make_tuple(id, ds->get_vol_card(id, bottom_level, top_level)));
	}
	return ret;
}
