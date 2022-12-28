//
//  Seq_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 24.03.22.
//

#include "Seq_Sketch.hpp"


//struct equation_data{
//	Eigen::SparseMatrix<double, Eigen::RowMajor> A_row;
//	double b;
//};
//
//
//double Seq_Sketch::f(unsigned int n, const double *x, double *grad, void *data) {
//
//	if (grad) {
//		for (int i = 0; i < n; i++) {
//			if (abs(x[i]) < 1e-10) {
//				grad[i] = 0.0;
//			} else if (x[i] < -1e-10) {
//				grad[i] = -1.0;
//			} else {
//				grad[i] = 1.0;
//			}
//		}
//	}
//	double ret = 0.0;
//	for (int i = 0; i < n; i ++) {
//		ret += abs(x[i]);
//	}
//	return ret;
//}
//
//double Seq_Sketch::eq(unsigned int n, const double *x, double *grad, void *data) {
//	Eigen::SparseMatrix<double, Eigen::RowMajor> A_row = ((equation_data*) data)->A_row;
//	double b = ((equation_data*) data)->b;
//	if (grad) {
//		memset((void*) grad, 0, n * sizeof(double));
//		for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_row, 0); it; ++it) {
//			double val = it.value();
//			int index = it.index();
//			grad[index] = val;
//		}
//	}
//	double ret = 0.0;
//	for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_row, 0); it; ++it) {
//		double val = it.value();
//		int index = it.index();
//		ret += val * x[index];
//	}
//	ret -= b;
//	return ret;
//}

struct func_data{
	Eigen::MatrixXd Nullsapce;
	Eigen::MatrixXd base;
};


double Seq_Sketch::f(unsigned int en, const double *x, double *grad, void *data) {
	
	func_data* func_data_ptr = (func_data*) data;
	
	Eigen::Map<Eigen::MatrixXd> z((double*) x, en, 1);
	
	Eigen::MatrixXd ret = func_data_ptr->base + func_data_ptr->Nullsapce * z;
	
	int m = (func_data_ptr->Nullsapce).rows();
	
	if (grad) {
		for (int i = 0; i < en; i++) {
			grad[i] = 0.0;
			for (int j = 0; j < m; j ++) {
				double multiplier = 0.0;
				if (ret(j, 0) < -1e-10) {
					multiplier = -1.0;
				} else if (ret(j, 0) > 1e-10) {
					multiplier = 1.0;
				}
				grad[i] += multiplier * (func_data_ptr->Nullsapce)(j, i);
			}
		}
	}
	
	return ret.array().abs().sum();
}

vector<double> Seq_Sketch::predict(arma::sp_dmat &hash_matrix, arma::dcolvec &volume_counters) {
	//	auto myfunc = [] (unsigned int n, const double* x, double* grad, void *data) {
	//		if (grad) {
	//			grad[0] = 2.0 * x[0];
	//			grad[1] = 2.0 * x[1];
	//		}
	//		return x[0] * x[0] + x[1] * x[1];
	//	};
	//
	//	double data[2][3] = {{1.0, 1.0, 2.0}, {1.0, -2.0, 2.0}};
	//
	//	auto myconstr = [] (unsigned int n, const double* x, double* grad, void *data) {
	//		if (grad) {
	//			grad[0] = ((double*) data)[0];
	//			grad[1] = ((double*) data)[1];
	//		}
	//		return x[0] * ((double*) data)[0] + x[1] * ((double*) data)[1] - ((double*) data)[2];
	//	};
	
	int n = (int) hash_matrix.n_cols;
	int m = (int) hash_matrix.n_rows;
	
//	arma::Mat<klab::DoubleReal> arma_hash_matrix(m, n, arma::fill::zeros);
	arma::Col<klab::DoubleReal> sol;
//	arma::Col<klab::DoubleReal> y(volume_counters.data(), volume_counters.size());
//	for (int k = 0; k < hash_matrix.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(hash_matrix, k); it; ++it) {
//			it.value();
//			it.row();   // row index
//			it.col();   // col index (here it is equal to k)
//			arma_hash_matrix(it.row(), it.col()) = it.value();
//		}
//	}
	arma::dmat arma_hash_matrix(hash_matrix);
	cout<<"Seq begin2"<<endl;
	
	klab::TSmartPointer<kl1p::TOperator<klab::DoubleReal> > op = new kl1p::TMatrixOperator<klab::DoubleReal>(arma_hash_matrix);
	kl1p::TBasisPursuitSolver<klab::DoubleReal> bp(1e-10);
	bp.solve(volume_counters, op, sol);
	
	vector<double>ret = vector<double>(sol.memptr(), sol.memptr() + n);
	//	cout<<hash_matrix<<endl;
	//	cout<<volume_counters<<endl;
	
	//	equation_data* data_array = (equation_data*) malloc(sizeof(equation_data) * m);
	//
	//	for (int i = 0; i < m; i++) {
	//		Eigen::SparseMatrix<double, Eigen::RowMajor> temp(hash_matrix.row(i));
	//		//		cout<<temp<<endl;
	//		new (&(data_array[i].A_row)) Eigen::SparseMatrix<double, Eigen::RowMajor>(hash_matrix.row(i));
	//		data_array[i].b = volume_counters(i, 0);
	//	}
	//
	//	double* lbs = (double*) calloc(n, sizeof(double));
	//
	//	auto PR_Sketch_f = [] (unsigned int n, const double *x, double *grad, void *data) {
	//		return f(n, x, grad, data);
	//	};
	//
	//	auto PR_Sketch_eq = [] (unsigned int n, const double *x, double *grad, void *data) {
	//		return eq(n, x, grad, data);
	//	};
	//
	//
	//
	//	nlopt_opt opt;
	//
	//	opt = nlopt_create(NLOPT_LD_SLSQP, n);
	////	opt = nlopt_create(NLOPT_GN_ISRES, n);
	//
	////	opt = nlopt_create(NLOPT_LN_COBYLA, n);
	//	nlopt_set_lower_bounds(opt, lbs);
	//	nlopt_set_min_objective(opt, PR_Sketch_f, NULL);
	//	for (int i = 0; i < m; i++) {
	//		if (data_array[i].A_row.nonZeros() > 0) {
	//			nlopt_add_equality_constraint(opt, PR_Sketch_eq, &(data_array[i]), 1e-10);
	//		}
	//	}
	//
	//	nlopt_set_xtol_rel(opt, 1e-10);
	//
	//	Eigen::MatrixXd densehashmatrix = hash_matrix.toDense();
	//	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomp(densehashmatrix);
	//	Eigen::MatrixXd intitial_guess = decomp.solve(volume_counters);
	//
	//	double* x = (double*) malloc(sizeof(double) * n);
	//	for (int i = 0; i < n; i++) {
	//		x[i] = intitial_guess(i, 0);
	//	}
	//	double minf;
	//
	//	if (nlopt_optimize(opt, x, &minf) < 0) {
	//		assert(false && "No solution found");
	//	}
	//
	//
	//
	//	vector<double> ret(x, x + n);
	//
	//	Eigen::Map<Eigen::MatrixXd> sol(ret.data(), n, 1);
	//	auto tempor = hash_matrix * intitial_guess;
	//	cout<<endl<<"Optimization error: "<<(volume_counters - tempor).norm()<<endl;
	//
	////	vector<double> ret(intitial_guess.data(), intitial_guess.data() + n);
	///
	///
	///
	///
	
	
	
	
	
	
	
	
	
//	Eigen::MatrixXd densehashmatrix = hash_matrix.toDense();
//	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
//	cod.compute(densehashmatrix);
//	Eigen::MatrixXd V = cod.matrixZ().transpose();
//	Eigen::MatrixXd Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
//	Eigen::MatrixXd P = cod.colsPermutation();
//	Null_space = P * Null_space; // Unpermute the columns
//
//	Eigen::MatrixXd initial_guess = cod.solve(volume_counters);
//
//	vector<double> ret;
//
//	if (Null_space.cols() == 0) {
//		auto tempor = hash_matrix * initial_guess;
//		cout<<endl<<"Optimization error: "<<(volume_counters - tempor).norm()<<endl;
//		ret = vector<double>(initial_guess.data(), initial_guess.data() + initial_guess.size());
//	} else {
//
//		auto Seq_Sketch_f = [] (unsigned int en, const double *x, double *grad, void *data) {
//			return f(en, x, grad, data);
//		};
//
//		func_data data;
//		data.base = initial_guess;
//		data.Nullsapce = Null_space;
//
//		nlopt_opt opt;
//
//		opt = nlopt_create(NLOPT_LD_SLSQP, Null_space.cols());
//		//	opt = nlopt_create(NLOPT_GN_ISRES, n);
//
////			opt = nlopt_create(NLOPT_LN_COBYLA, Null_space.cols());
//		nlopt_set_min_objective(opt, Seq_Sketch_f, (void*) &data);
//
//		nlopt_set_xtol_rel(opt, 1e-10);
//
//		double* z = (double*) malloc(sizeof(double) * Null_space.cols());
//		for (int i = 0; i < Null_space.cols(); i++) {
//			z[i] = 0;
//		}
//		double minf;
//
//		if (nlopt_optimize(opt, z, &minf) < 0) {
//			assert(false && "No solution found");
//		}
//
//		Eigen::Map<Eigen::MatrixXd> zmat((double*) z, Null_space.cols(), 1);
//
//		Eigen::MatrixXd sol = initial_guess + Null_space * zmat;
//
//		auto tempor = hash_matrix * sol;
//		cout<<endl<<"Optimization error: "<<(volume_counters - tempor).norm()<<endl;
//
//
//		ret = vector<double>(sol.data(), sol.data() + sol.size());
//
//	}
	
	return ret;
}
