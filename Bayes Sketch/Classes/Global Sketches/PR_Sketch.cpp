//
//  PR_Sketch.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 23.03.22.
//

#include "PR_Sketch.hpp"

struct equation_data{
	Eigen::SparseMatrix<double, Eigen::RowMajor> A_row;
	double b;
};


double PR_Sketch::f(unsigned int n, const double *x, double *grad, void *data) {
	if (grad) {
		memcpy((void*) grad, (const void*) x, n * sizeof(double));
	}
	double ret = 0.0;
	for (int i = 0; i < n; i ++) {
		ret += x[i] * x[i];
	}
	return ret * 0.5;
}

double PR_Sketch::eq(unsigned int n, const double *x, double *grad, void *data) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> A_row = ((equation_data*) data)->A_row;
	double b = ((equation_data*) data)->b;
	if (grad) {
		memset((void*) grad, 0, n * sizeof(double));
		for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_row, 0); it; ++it) {
			double val = it.value();
			int index = it.index();
			grad[index] = val;
		}
	}
	double ret = 0.0;
	for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_row, 0); it; ++it) {
		double val = it.value();
		int index = it.index();
		ret += val * x[index];
	}
	ret -= b;
	return ret;
}

vector<double> PR_Sketch::predict(arma::sp_dmat &hash_matrix, arma::dcolvec &volume_counters) {
////	auto myfunc = [] (unsigned int n, const double* x, double* grad, void *data) {
////		if (grad) {
////			grad[0] = 2.0 * x[0];
////			grad[1] = 2.0 * x[1];
////		}
////		return x[0] * x[0] + x[1] * x[1];
////	};
////
////	double data[2][3] = {{1.0, 1.0, 2.0}, {1.0, -2.0, 2.0}};
////
////	auto myconstr = [] (unsigned int n, const double* x, double* grad, void *data) {
////		if (grad) {
////			grad[0] = ((double*) data)[0];
////			grad[1] = ((double*) data)[1];
////		}
////		return x[0] * ((double*) data)[0] + x[1] * ((double*) data)[1] - ((double*) data)[2];
////	};
//
//	int n = (int) hash_matrix.cols();
//	int m = (int) hash_matrix.rows();
//
////	cout<<hash_matrix.transpose()<<endl;
//
//	equation_data* data_array = (equation_data*) malloc(sizeof(equation_data) * m);
//
//	for (int i = 0; i < m; i++) {
//		Eigen::SparseMatrix<double, Eigen::RowMajor> temp(hash_matrix.row(i));
////		cout<<temp<<endl;
//		auto t = new (&(data_array[i].A_row)) Eigen::SparseMatrix<double, Eigen::RowMajor>(hash_matrix.row(i));
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
//	nlopt_set_lower_bounds(opt, lbs);
//	nlopt_set_min_objective(opt, PR_Sketch_f, NULL);
//	for (int i = 0; i < m; i++) {
//		nlopt_add_equality_constraint(opt, PR_Sketch_eq, &(data_array[i]), 1e-10);
//	}
//
//	nlopt_set_xtol_rel(opt, 1e-10);
//	double* x = (double*) malloc(sizeof(double) * n);
//	double init_est = l1 / n;
//	for (int i = 0; i < n; i++) {
//		x[i] = init_est;
//	}
//	double minf;
//	if (nlopt_optimize(opt, x, &minf) < 0) {
//		printf("nlopt failed!\n");
//	}
//	else {
//		printf("found minimum \n");
//	}
//
//
//
//	vector<double> ret(x, x + n);
//
//	Eigen::Map<Eigen::MatrixXd> sol(ret.data(), n, 1);
	
	
	
//	 Skeleton for Sparse Pseudoinverse comoutation. Not finished because it takes longer than dense QR.
//	Eigen::SparseMatrix<double, Eigen::ColMajor> colmajor_hash_matrix(hash_matrix);
//	colmajor_hash_matrix.makeCompressed();
//
//	Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::COLAMDOrdering<int>> qr(colmajor_hash_matrix);
//	Eigen::MatrixXd qr_sol = qr.solve(volume_counters);

//	auto Qmat = qr.matrixQ();
//	auto Rmat = qr.matrixR().triangularView<Eigen::Upper>();
//	Eigen::MatrixXd intermediate = Qmat.transpose() * volume_counters;
//	cout<<Rmat.cols()<<endl<<Rmat.rows()<<endl;
//	Eigen::MatrixXd almostsol = Rmat.solve(intermediate);
//	Eigen::MatrixXd qr_sol = qr.colsPermutation() * almostsol;

	
//	auto temporqr = hash_matrix * qr_sol;
//	cout<<endl<<"PR Sketch optimization error: "<<(volume_counters - temporqr).norm()<<endl;
	
	
	int n = (int) hash_matrix.n_cols;
	int m = (int) hash_matrix.n_rows;
	
//	arma::umat locs(2, hash_matrix.nonZeros());
//	arma::dcolvec vals(hash_matrix.nonZeros());
//	int ind = 0;
//	for (int k = 0; k < hash_matrix.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(hash_matrix, k); it; ++it) {
//			it.value();
//			it.row();   // row index
//			it.col();   // col index (here it is equal to k)
//			locs(0, ind) = it.row();
//			locs(1, ind) = it.col();
//			vals(ind) = it.value();
//			ind++;
//		}
//	}
//
//	arma::sp_dmat arma_hash_matrix(locs, vals);
//	arma::dcolvec arma_volume_counters(volume_counters.data(), volume_counters.size());
//
	arma::dmat U;
	arma::dcolvec s;
	arma::dmat V;
	
	arma::svds(U, s, V, hash_matrix, volume_counters.n_elem);
	
	int rank = 0;
	for (; rank < volume_counters.n_elem; rank++){
		if (s(rank) <= 1e-10) {
			break;
		}
	}
//	cout<<U.n_rows<<" "<<U.n_cols<<endl;
//	cout<<V.n_rows<<" "<<V.n_cols<<endl;
//	cout<<rank<<endl;
//	cout<<m<<endl;
//	cout<<n<<endl;
	U = U.submat(0, 0, m - 1, rank - 1);
	V = V.submat(0, 0, n - 1, rank - 1);
	s = s.head(rank);
	
	
	
	arma::dmat Sinv =  arma::diagmat(1.0 / s);
//	cout<<U<<endl<<s<<endl<<V<<endl<<arma_hash_matrix<<endl<<Sinv<<endl;
//	cout<<s<<endl<<(1.0 / s)<<endl;
//
//	cout<<Sinv<<endl;
	
//	cout<<U.n_rows<<" "<<U.n_cols<<endl;
//		cout<<V.n_rows<<" "<<V.n_cols<<endl;
//	cout<<Sinv.n_rows<<" "<<Sinv.n_cols<<endl;
	
	arma::dcolvec sol = V * Sinv * U.t() * volume_counters;
//	cout<<sol.size()<<endl;
//	cout<<sol<<endl;
	
	
	

	
	
	
	
	
	
	
	
	


//	// Compute pseudoinverse from dense hash matrix. Is not memory efficient at all but it appears to be the fastest.
//	Eigen::MatrixXd densehashmatrix = hash_matrix.toDense();
//	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomp(densehashmatrix);
//	Eigen::MatrixXd sol = decomp.solve(volume_counters);
////
////	cout<<sol<<endl;
////
//
//	auto tempor = hash_matrix * sol;
//	cout<<endl<<"PR Sketch optimization error: "<<(volume_counters - tempor).norm()<<endl;
//
////	cout<<(sol-qr_sol).norm()<<endl;
	
//	vector<double> ret(sol.data(), sol.data() + sol.size());
	vector<double> ret(sol.memptr(), sol.memptr() + n);
	return ret;
}
