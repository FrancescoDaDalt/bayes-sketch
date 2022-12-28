//
//  main.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 16.03.22.
//

#include <iostream>


#include "Classes/Local Sketches/C_Sketch.hpp"
#include "Classes/Local Sketches/CM_Sketch.hpp"
#include "Classes/Local Sketches/CB_Sketch.hpp"
#include "Classes/Local Sketches/CCB_Sketch.hpp"
#include "Classes/Local Sketches/CCA_Sketch.hpp"
#include "Classes/Global Sketches/PR_Sketch.hpp"
#include "Classes/Global Sketches/Seq_Sketch.hpp"
#include "Classes/Trace_Container.hpp"
#include "Classes/Datastructure Interfaces/LocalDatastructure_Interface.hpp"
#include "Classes/Sketch Interfaces/LocalSketch_Interface.hpp"
#include "Classes/Sketch Interfaces/GlobalSketch_Interface.hpp"
#include "Classes/Sketch Utility/LocalSketch_Engine.hpp"
#include "Classes/Sketch Utility/GlobalSketch_Engine.hpp"
#include "Classes/Datastructures/SingleHash.hpp"
#include "Classes/Datastructures/DoubleHash_Centered.hpp"
#include "Classes/Datastructures/DoubleHash_Normalized.hpp"
#include "Classes/Dump_Container.hpp"
//#include "CompressedSensingExample.h"


int main(int argc, const char * argv[]) {
	
	// Construct Containers that will dump the statistics into a csv file
	Dump_Container c_dump("C Sketch Statistics", 20);
	Dump_Container cm_dump("CM Sketch Statistics", 20);
//	Dump_Container cb_dump("CB Sketch Statistics", 20);
//	Dump_Container ccb_dump("CCB Sketch Statistics", 20);
	Dump_Container cbtrainedperf_dump("Trained rsquare_surrogate perfect mu CB Sketch Statistics", 20);
	Dump_Container ccbtrainedperf_dump("Trained rsquare_surrogate perfect mu CCB Sketch Statistics", 20);
//	Dump_Container cbtrainedzero_dump("Trained avg_rel_error zero mu CB Sketch Statistics", 20);
//	Dump_Container ccbtrainedzero_dump("Trained avg_rel_error zero mu CCB Sketch Statistics", 20);
//	Dump_Container cca_dump("CCA Sketch Statistics", 20);
	Dump_Container seq_dump("Seq Sketch Statistics", 20);
	Dump_Container pr_dump("PR Sketch Statistics", 20);
	
	// Set random seed and get random seed for trace synthesizaion
	random_device rd;
	auto seed = rd();
	srand(time(NULL));
	
	// Read in trace froma csv file with the following columns : (timestamp, key, payload).
	// The timestamp is ignored and the data stream is assumed to be in the right ordering.
	// This constructor subsamles randomly 10000 keys from the stream. The random seed is given by seed.
//	Trace_Container trace("Data/448000-1.csv", 10000, seed);
//	Trace_Container train_trace("Data/448000-1.csv", 10000, seed + 42);
	
	// This constructor subsamples the data stream by extracting the segment specified by the tuple.
	// In this example, the trainign trace takes the first 10% and the test trace the last 10% of teh datastream.
	// num_bloom_filter_cells specifies how many bloom filter cells to use for cardinality estimation. Whether the filter is used or not can be specified later.
	int num_bloom_filter_cells = 100000;
	Trace_Container trace("Data/univ/univ1_5.csv", tuple<float, float>(0.9, 1.0), num_bloom_filter_cells);
	Trace_Container train_trace("Data/univ/univ1_1.csv", tuple<float, float>(0.0, 0.1), num_bloom_filter_cells);
	
	// Sythesize trace based on a distribution
//	Trace_Container trace(10000);
//	trace.synthesizeTranslatedPoisson(100, 100, seed);
//	Trace_Container train_trace(10000);
//	train_trace.synthesizeTranslatedPoisson(100, 100, seed + 42);
//	Trace_Container trace(10000);
//	trace.synthesizeTranslatedGeometric(100, 100, seed);
//	Trace_Container train_trace(10000);
//	train_trace.synthesizeTranslatedGeometric(100, 100, seed + 42);
//	Trace_Container trace(100000);
//	trace.synthesizeTranslatedStudentT(1, 2000, 100, seed);
	
	// Define the "grid" of datastructure dimensions we want to try
	// In this example, we only want to do one run with a 30x30 datastructure.
	vector<int> wvals{30};
	vector<int> dvals{30};
	
	// Specify the artificial noise level to be used when estimating the cardinalitis of counters.
	// 1.0 is maximum and 0.0 is no noise.
	double noise_level = 0.0;
	
	// Whether to use the bloom filter when computing cardinality estimates. The trace that is passed to the algorithm must have been created with a bloom filter for this to work.
	// Artificial noise and bloom-filter noise cannot be simulated at the same time.
	bool use_bloom_filter = true;
	for (int w: wvals) {
		for (int d : dvals) {
			
			// Feed traces to all datastructures used. Datastructures for the PR- and Seq-Sketch have depth 5 as they tend to perform better in these circumstances.
			SingleHash singlehash_ds(d, w);
			singlehash_ds.feed_trace(&trace, use_bloom_filter, noise_level);

			SingleHash singlehash_train_ds(d, w);
			singlehash_train_ds.feed_trace(&train_trace, use_bloom_filter, noise_level);
//
			SingleHash singlehash_pr_ds(5, d * w / 5);
			singlehash_pr_ds.feed_trace(&trace, use_bloom_filter, noise_level);
//
			DoubleHash_Centered doublehash_ds(d, w);
			doublehash_ds.feed_trace(&trace, use_bloom_filter, noise_level);
			
			DoubleHash_Normalized doublehash_norm_ds(5, d * w / 5);
			doublehash_norm_ds.feed_trace(&trace, use_bloom_filter, noise_level);
			
			// Construct all the sketches
			C_Sketch c_sketch;
			CM_Sketch cm_sketch;
//			CB_Sketch cb_sketch;
//			CCB_Sketch ccb_sketch;
			CB_Sketch cbtrainedperf_sketch;
			CCB_Sketch ccbtrainedperf_sketch;
//			CB_Sketch cbtrainedzero_sketch;
//			CCB_Sketch ccbtrainedzero_sketch;
//			CCA_Sketch cca_sketch;
			Seq_Sketch seq_sketch;
			PR_Sketch pr_sketch;
			
			// Define all prior chi parameters to try out during training
			// This is just one of uncountable ways how optimal priors may be found. Finding the optimal prior by menas of bayesian optimization is likely faster but a bit of an overkill.
			// As described in the paper, we use a log-scaled grid to find the best chi. The search range comes from an educated guess of where the optimal chi should lie based on the
			// definition of chi.
			double basechi = 1.0 / ((singlehash_train_ds.get_l1() * (double) singlehash_train_ds.get_l1()) / train_trace.get_distinct_ids().size());
			vector<double> chivalues;
			for (int i = 0; i < 200; i++) {
				chivalues.push_back(basechi * pow(1.12, (double) i));
//				chivalues.push_back(basechi * pow(1.254, (double) i));
			}
//			// Define all prior mu parameters to try out
			vector<double> muvalues_perf{((double) singlehash_train_ds.get_l1()) / train_trace.get_distinct_ids().size()};
			vector<double> muvalues_zero{0.0};
			auto trainids = train_trace.get_distinct_ids();
			auto traingt = train_trace.get_ground_truth();
//
//			// Perform "training". rsquare_surrogate is the optimization goal
			cbtrainedperf_sketch.minimize_on_trace((LocalDatastructure_Interface*) &singlehash_train_ds, trainids, traingt, chivalues, muvalues_perf, "rsquare_surrogate");
			ccbtrainedperf_sketch.minimize_on_trace((LocalDatastructure_Interface*) &singlehash_train_ds, trainids, traingt, chivalues, muvalues_perf, "rsquare_surrogate");
//			cbtrainedzero_sketch.minimize_on_trace((LocalDatastructure_Interface*) &singlehash_train_ds, trainids, traingt, chivalues, muvalues_zero, "avg_rel_error");
//			ccbtrainedzero_sketch.minimize_on_trace((LocalDatastructure_Interface*) &singlehash_train_ds, trainids, traingt, chivalues, muvalues_zero, "avg_rel_error");
//
//			// Initialize all the metadata that we want to dump together with the error statistics
			vector<tuple<string, string>> c_met;
			c_met.push_back(make_tuple("Algorithm", "Count Sketch"));
			c_met.push_back(make_tuple("Width", to_string(w)));
			c_met.push_back(make_tuple("Depth", to_string(d)));
			c_met.push_back(make_tuple("Memory", to_string(d * w)));
			c_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			c_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			c_met.push_back(make_tuple("Trace", trace.get_trace()));

			vector<tuple<string, string>> cm_met;
			cm_met.push_back(make_tuple("Algorithm", "CM Sketch"));
			cm_met.push_back(make_tuple("Width", to_string(w)));
			cm_met.push_back(make_tuple("Depth", to_string(d)));
			cm_met.push_back(make_tuple("Memory", to_string(d * w)));
			cm_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			cm_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			cm_met.push_back(make_tuple("Trace", trace.get_trace()));
//
//			vector<tuple<string, string>> cb_met;
//			cb_met.push_back(make_tuple("Algorithm", "CB Sketch"));
//			cb_met.push_back(make_tuple("Prior", cb_sketch.get_prior()));
//			cb_met.push_back(make_tuple("Mu", cb_sketch.get_priormu()));
//			cb_met.push_back(make_tuple("Chi", cb_sketch.get_priorchi()));
//			cb_met.push_back(make_tuple("Width", to_string(w)));
//			cb_met.push_back(make_tuple("Depth", to_string(d)));
//			cb_met.push_back(make_tuple("Memory", to_string(d * w)));
//			cb_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
//			cb_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
//			cb_met.push_back(make_tuple("Trace", trace.get_trace()));
//
//			vector<tuple<string, string>> ccb_met;
//			ccb_met.push_back(make_tuple("Algorithm", "CCB Sketch"));
//			ccb_met.push_back(make_tuple("Prior", ccb_sketch.get_prior()));
//			ccb_met.push_back(make_tuple("Mu", ccb_sketch.get_priormu()));
//			ccb_met.push_back(make_tuple("Chi", ccb_sketch.get_priorchi()));
//			ccb_met.push_back(make_tuple("Width", to_string(w)));
//			ccb_met.push_back(make_tuple("Depth", to_string(d)));
//			ccb_met.push_back(make_tuple("Memory", to_string(d * w)));
//			ccb_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
//			ccb_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
//			ccb_met.push_back(make_tuple("Trace", trace.get_trace()));
//			ccb_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
//
			vector<tuple<string, string>> cbtrainedperf_met;
			cbtrainedperf_met.push_back(make_tuple("Algorithm", "CB Sketch"));
			cbtrainedperf_met.push_back(make_tuple("Prior", cbtrainedperf_sketch.get_prior()));
			cbtrainedperf_met.push_back(make_tuple("Mu", cbtrainedperf_sketch.get_priormu()));
			cbtrainedperf_met.push_back(make_tuple("Chi", cbtrainedperf_sketch.get_priorchi()));
			cbtrainedperf_met.push_back(make_tuple("Width", to_string(w)));
			cbtrainedperf_met.push_back(make_tuple("Depth", to_string(d)));
			cbtrainedperf_met.push_back(make_tuple("Memory", to_string(d * w)));
			cbtrainedperf_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			cbtrainedperf_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			cbtrainedperf_met.push_back(make_tuple("Trace", trace.get_trace()));

			vector<tuple<string, string>> ccbtrainedperf_met;
			ccbtrainedperf_met.push_back(make_tuple("Algorithm", "CCB Sketch"));
			ccbtrainedperf_met.push_back(make_tuple("Prior", ccbtrainedperf_sketch.get_prior()));
			ccbtrainedperf_met.push_back(make_tuple("Mu", ccbtrainedperf_sketch.get_priormu()));
			ccbtrainedperf_met.push_back(make_tuple("Chi", ccbtrainedperf_sketch.get_priorchi()));
			ccbtrainedperf_met.push_back(make_tuple("Width", to_string(w)));
			ccbtrainedperf_met.push_back(make_tuple("Depth", to_string(d)));
			ccbtrainedperf_met.push_back(make_tuple("Memory", to_string(d * w)));
			ccbtrainedperf_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			ccbtrainedperf_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			ccbtrainedperf_met.push_back(make_tuple("Trace", trace.get_trace()));
			ccbtrainedperf_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
			ccbtrainedperf_met.push_back(make_tuple("Bloom Filter", to_string(use_bloom_filter)));
			ccbtrainedperf_met.push_back(make_tuple("Bloom Filter Size", to_string(num_bloom_filter_cells)));
//
//			vector<tuple<string, string>> cbtrainedzero_met;
//			cbtrainedzero_met.push_back(make_tuple("Algorithm", "CB Sketch"));
//			cbtrainedzero_met.push_back(make_tuple("Prior", cbtrainedzero_sketch.get_prior()));
//			cbtrainedzero_met.push_back(make_tuple("Mu", cbtrainedzero_sketch.get_priormu()));
//			cbtrainedzero_met.push_back(make_tuple("Chi", cbtrainedzero_sketch.get_priorchi()));
//			cbtrainedzero_met.push_back(make_tuple("Width", to_string(w)));
//			cbtrainedzero_met.push_back(make_tuple("Depth", to_string(d)));
//			cbtrainedzero_met.push_back(make_tuple("Memory", to_string(d * w)));
//			cbtrainedzero_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
//			cbtrainedzero_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
//			cbtrainedzero_met.push_back(make_tuple("Trace", trace.get_trace()));
//
//			vector<tuple<string, string>> ccbtrainedzero_met;
//			ccbtrainedzero_met.push_back(make_tuple("Algorithm", "CCB Sketch"));
//			ccbtrainedzero_met.push_back(make_tuple("Prior", ccbtrainedzero_sketch.get_prior()));
//			ccbtrainedzero_met.push_back(make_tuple("Mu", ccbtrainedzero_sketch.get_priormu()));
//			ccbtrainedzero_met.push_back(make_tuple("Chi", ccbtrainedzero_sketch.get_priorchi()));
//			ccbtrainedzero_met.push_back(make_tuple("Width", to_string(w)));
//			ccbtrainedzero_met.push_back(make_tuple("Depth", to_string(d)));
//			ccbtrainedzero_met.push_back(make_tuple("Memory", to_string(d * w)));
//			ccbtrainedzero_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
//			ccbtrainedzero_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
//			ccbtrainedzero_met.push_back(make_tuple("Trace", trace.get_trace()));
//			ccbtrainedzero_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
//
//			vector<tuple<string, string>> cca_met;
//			cca_met.push_back(make_tuple("Algorithm", "CCA Sketch"));
//			cca_met.push_back(make_tuple("Width", to_string(w)));
//			cca_met.push_back(make_tuple("Depth", to_string(d)));
//			cca_met.push_back(make_tuple("Memory", to_string(d * w)));
//			cca_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
//			cca_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
//			cca_met.push_back(make_tuple("Trace", trace.get_trace()));
//			cca_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
			
			vector<tuple<string, string>> seq_met;
			seq_met.push_back(make_tuple("Algorithm", "Seq Sketch"));
			seq_met.push_back(make_tuple("Width", to_string(d * w / 5)));
			seq_met.push_back(make_tuple("Depth", to_string(5)));
			seq_met.push_back(make_tuple("Memory", to_string(d * w)));
			seq_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			seq_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			seq_met.push_back(make_tuple("Trace", trace.get_trace()));
			seq_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
			seq_met.push_back(make_tuple("Bloom Filter", to_string(use_bloom_filter)));
			seq_met.push_back(make_tuple("Bloom Filter Size", to_string(num_bloom_filter_cells)));
			
			vector<tuple<string, string>> pr_met;
			pr_met.push_back(make_tuple("Algorithm", "PR Sketch"));
			pr_met.push_back(make_tuple("Width", to_string(w * d / 5)));
			pr_met.push_back(make_tuple("Depth", to_string(5)));
			pr_met.push_back(make_tuple("Memory", to_string(d * w)));
			pr_met.push_back(make_tuple("Trace Type", trace.get_tracetype()));
			pr_met.push_back(make_tuple("Num Items", to_string(trace.get_distinct_ids().size())));
			pr_met.push_back(make_tuple("Trace", trace.get_trace()));
			pr_met.push_back(make_tuple("Cardinality Noise", to_string(noise_level)));
			pr_met.push_back(make_tuple("Bloom Filter", to_string(use_bloom_filter)));
			pr_met.push_back(make_tuple("Bloom Filter Size", to_string(num_bloom_filter_cells)));
//
//			// For all sketches, create an execution engine based on datastructure and query algorithm
//			// Then, let the sketch engine compute the estimate of the aggregate key volumes
//			// And insert the computed error statistics into the dump container
			LocalSketch_Engine c_sketch_engine((LocalSketch_Interface*) &c_sketch,(LocalDatastructure_Interface*) &doublehash_ds);
			vector<tuple<string, double>>  c_results = c_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 500000);
			c_dump.insert_row(c_results, c_met);

			LocalSketch_Engine cm_sketch_engine((LocalSketch_Interface*) &cm_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
			vector<tuple<string, double>>  cm_results = cm_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 500000);
			cm_dump.insert_row(cm_results, cm_met);
//
//			LocalSketch_Engine cb_sketch_engine((LocalSketch_Interface*) &cb_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
//			vector<tuple<string, double>>  cb_results = cb_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 1000000);
//			cb_dump.insert_row(cb_results, cb_met);
//
			LocalSketch_Engine cbtrainedperf_sketch_engine((LocalSketch_Interface*) &cbtrainedperf_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
			vector<tuple<string, double>>  cbtrainedperf_results = cbtrainedperf_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 500000);
			cbtrainedperf_dump.insert_row(cbtrainedperf_results, cbtrainedperf_met);

//			LocalSketch_Engine cbtrainedzero_sketch_engine((LocalSketch_Interface*) &cbtrainedzero_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
//			vector<tuple<string, double>>  cbtrainedzero_results = cbtrainedzero_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 1000000);
//			cbtrainedzero_dump.insert_row(cbtrainedzero_results, cbtrainedzero_met);
//
//			LocalSketch_Engine ccb_sketch_engine((LocalSketch_Interface*) &ccb_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
//			vector<tuple<string, double>>  ccb_results = ccb_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 1000000);
//			ccb_dump.insert_row(ccb_results, ccb_met);
//
			LocalSketch_Engine ccbtrainedperf_sketch_engine((LocalSketch_Interface*) &ccbtrainedperf_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
			vector<tuple<string, double>>  ccbtrainedperf_results = ccbtrainedperf_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 500000);
			ccbtrainedperf_dump.insert_row(ccbtrainedperf_results, ccbtrainedperf_met);

//			LocalSketch_Engine ccbtrainedzero_sketch_engine((LocalSketch_Interface*) &ccbtrainedzero_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
//			vector<tuple<string, double>>  ccbtrainedzero_results = ccbtrainedzero_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 1000000);
//			ccbtrainedzero_dump.insert_row(ccbtrainedzero_results, ccbtrainedzero_met);
//
//			LocalSketch_Engine cca_sketch_engine((LocalSketch_Interface*) &cca_sketch,(LocalDatastructure_Interface*) &singlehash_ds);
//			vector<tuple<string, double>>  cca_results = cca_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth(), 8, 1000000);
//			cca_dump.insert_row(cca_results, cca_met);
			
//			cout<<"Seq begin"<<endl;
			GlobalSketch_Engine seq_sketch_engine((GlobalSketch_Interface*) &seq_sketch,(GlobalDatastructure_Interface*) &doublehash_norm_ds);
			vector<tuple<string, double>>  seq_results = seq_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth());
			seq_dump.insert_row(seq_results, seq_met);
//			cout<<"Seq end"<<endl;

//			cout<<"PR begin"<<endl;
			GlobalSketch_Engine pr_sketch_engine((GlobalSketch_Interface*) &pr_sketch,(GlobalDatastructure_Interface*) &singlehash_pr_ds);
			vector<tuple<string, double>>  pr_results = pr_sketch_engine.predict_and_benchmark(trace.get_distinct_ids(), trace.get_ground_truth());
			pr_dump.insert_row(pr_results, pr_met);
//			cout<<"PR end"<<endl;
//
//			// Print some statistics of interest
			cout<<endl<<"C Sketch:"<<endl;
			for (tuple<string, double> result : c_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
			cout<<endl<<"CM Sketch:"<<endl;
			for (tuple<string, double> result : cm_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
//			cout<<endl<<"CB Sketch:"<<endl;
//			for (tuple<string, double> result : cb_results) {
//				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
//					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
//				}
//			}
			cout<<endl<<"Trained Perf CB Sketch:"<<endl;
			for (tuple<string, double> result : cbtrainedperf_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0 ) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
//			cout<<endl<<"Trained Zero CB Sketch:"<<endl;
//			for (tuple<string, double> result : cbtrainedzero_results) {
//				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
//					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
//				}
//			}
//			cout<<endl<<"CCA Sketch:"<<endl;
//			for (tuple<string, double> result : cca_results) {
//				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
//					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
//				}
//			}
//			cout<<endl<<"CCB Sketch:"<<endl;
//			for (tuple<string, double> result : ccb_results) {
//				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
//					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
//				}
//			}
			cout<<endl<<"Trained Perf CCB Sketch:"<<endl;
			for (tuple<string, double> result : ccbtrainedperf_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
//			cout<<endl<<"Trained Zero CCB Sketch:"<<endl;
//			for (tuple<string, double> result : ccbtrainedzero_results) {
//				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
//					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
//				}
//			}
			
			
			cout<<endl<<"Seq Sketch:"<<endl;
			for (tuple<string, double> result : seq_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
			cout<<endl<<"PR Sketch:"<<endl;
			for (tuple<string, double> result : pr_results) {
				if (get<0>(result).compare("avg_rel_error") == 0 || get<0>(result).compare("avg_abs_error") == 0 || get<0>(result).compare("rsquare_surrogate") == 0) {
					cout<<"Error statistic: "<<get<0>(result)<<" \t Value: "<<get<1>(result)<<endl;
				}
			}
			
		}
	}
	
	// Dump the simulation results as csv file
	c_dump.dump_to_csv();
	cm_dump.dump_to_csv();
//	cb_dump.dump_to_csv();
//	ccb_dump.dump_to_csv();
	cbtrainedperf_dump.dump_to_csv();
	ccbtrainedperf_dump.dump_to_csv();
//	cbtrainedzero_dump.dump_to_csv();
//	ccbtrainedzero_dump.dump_to_csv();
//	cca_dump.dump_to_csv();
	seq_dump.dump_to_csv();
	pr_dump.dump_to_csv();
	
	return 0;
}
