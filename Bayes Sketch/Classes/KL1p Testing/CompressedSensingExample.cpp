// KL1p - A portable C++ compressed sensing library.
// Copyright (c) 2011-2012 René Gebel
// 
// This file is part of the KL1p C++ library.
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY of fitness for any purpose. 
//
// This library is free software; You can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License (LGPL) 
// as published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// See http://www.opensource.org/licenses for more info.

#include "CompressedSensingExample.h"

using namespace kl1p;




// ---------------------------------------------------------------------------------------------------- //

void	kl1p::CreateGaussianSignal(klab::UInt32 size, klab::UInt32 sparsity, klab::DoubleReal mean, klab::DoubleReal sigma, arma::Col<klab::DoubleReal>& out)
{
	out.set_size(size);
	out.fill(0.0);

	std::vector<klab::TArrayElement<klab::DoubleReal> > indices;
    for(klab::UInt32 i=0; i<size; ++i)
        indices.push_back(klab::TArrayElement<klab::DoubleReal>(i, klab::KRandom::Instance().generateDoubleReal(0.0, 1.0)));

    std::partial_sort(indices.begin(), indices.begin()+klab::Min(size, sparsity), indices.end(), std::greater<klab::TArrayElement<klab::DoubleReal> >());  

	for(klab::UInt32 i=0; i<sparsity; ++i)
	{
		klab::DoubleReal u1 = klab::KRandom::Instance().generateDoubleReal(0.0, 1.0);
		klab::DoubleReal u2 = klab::KRandom::Instance().generateDoubleReal(0.0, 1.0);

		klab::DoubleReal sign = klab::KRandom::Instance().generateBool() ? -1.0 : 1.0;
		out[indices[i].i()] = sign * ((klab::Sqrt(-2.0*klab::Log(u1)) * klab::Cos(2.0*klab::PI*u2))*sigma + mean);
	}
}

// ---------------------------------------------------------------------------------------------------- //

void	kl1p::WriteToCSVFile(const arma::Col<klab::DoubleReal>& signal, const std::string& filePath)
{
	std::ofstream of(filePath.c_str());
	if(of.is_open())
	{
		for(klab::UInt32 i=0; i<signal.n_rows; ++i)
			of<<i<<";"<<signal[i]<<std::endl;

		of.close();
	}
	else
	{
		std::cout<<"ERROR! Unable to open file \""<<filePath<<"\" !"<<std::endl;
	}
}

// ---------------------------------------------------------------------------------------------------- //

void	kl1p::RunExample()
{
	try
	{
		
		
		arma::Mat<klab::DoubleReal> test_matrix(2, 2);
		arma::Col<klab::DoubleReal> x0(2);
		arma::Col<klab::DoubleReal> x;
		arma::Col<klab::DoubleReal> y;
		test_matrix(0, 0) = 1.0;
		test_matrix(0, 1) = 2.0;
		test_matrix(1, 0) = 3.0;
		test_matrix(1, 1) = 4.0;
		
		x0(0) = 1.0;
		x0(1) = 2.0;
		klab::TSmartPointer<kl1p::TOperator<klab::DoubleReal> > op = new kl1p::TMatrixOperator<klab::DoubleReal>(test_matrix);
		op->apply(x0, y);
		std::cout<<y(0)<<std::endl<<y(1)<<std::endl;
		kl1p::TBasisPursuitSolver<klab::DoubleReal> bp(0.01);
		bp.solve(y, op, x);
		
		std::cout<<x0<<std::endl;
		std::cout<<x<<std::endl;
		
		
//
//
//
//
//		std::cout<<"Start of KL1p compressed-sensing example."<<std::endl;
//		std::cout<<"Try to determine a sparse vector x "<<std::endl;
//		std::cout<<"from an underdetermined set of linear measurements y=A*x, "<<std::endl;
//		std::cout<<"where A is a random gaussian i.i.d sensing matrix."<<std::endl;
//
//		klab::UInt32 n = 256;					// Size of the original signal x0.
//		klab::DoubleReal alpha = 0.5;			// Ratio of the cs-measurements.
//		klab::DoubleReal rho = 0.1;				// Ratio of the sparsity of the signal x0.
//		klab::UInt32 m = klab::UInt32(alpha*n);	// Number of cs-measurements.
//		klab::UInt32 k = klab::UInt32(rho*n);	// Sparsity of the signal x0 (number of non-zero elements).
//		klab::UInt64 seed = 0;					// Seed used for random number generation (0 if regenerate random numbers on each launch).
//		bool bWrite = false;					// Write signals to files ?
//
//		// Initialize random seed if needed.
//		if(seed > 0)
//			klab::KRandom::Instance().setSeed(seed);
//
//		// Display signal informations.
//		std::cout<<"=============================="<<std::endl;
//		std::cout<<"N="<<n<<" (signal size)"<<std::endl;
//		std::cout<<"M="<<m<<"="<<std::setprecision(5)<<(alpha*100.0)<<"% (number of measurements)"<<std::endl;
//		std::cout<<"K="<<k<<"="<<std::setprecision(5)<<(rho*100.0)<<"% (signal sparsity)"<<std::endl;
//		std::cout<<"Random Seed="<<klab::KRandom::Instance().seed()<<std::endl;
//		std::cout<<"=============================="<<std::endl;
//
//		arma::Col<klab::DoubleReal> x0;					// Original signal x0 of size n.
//		kl1p::CreateGaussianSignal(n, k, 0.0, 1.0, x0);	// Create randomly the original signal x0.
//
//		if(bWrite)
//			kl1p::WriteToCSVFile(x0, "OriginalSignal.csv");	// Write x0 to a file.
//
//		// Create random gaussian i.i.d matrix A of size (m,n).
//		klab::TSmartPointer<kl1p::TOperator<klab::DoubleReal> > A = new kl1p::TNormalRandomMatrixOperator<klab::DoubleReal>(m, n, 0.0, 1.0);
//		A  = new kl1p::TScalingOperator<klab::DoubleReal>(A, 1.0/klab::Sqrt(klab::DoubleReal(m)));	// Pseudo-normalization of the matrix (required for AMP and EMBP solvers).
//
//		// Perform cs-measurements of size m.
//		arma::Col<klab::DoubleReal> y;
//		A->apply(x0, y);
//
//		klab::DoubleReal tolerance = 1e-3;	// Tolerance of the solution.
//		arma::Col<klab::DoubleReal> x;		// Will contain the solution of the reconstruction.
//
//		klab::KTimer timer;
//
//		// Compute Basis-Pursuit.
//		std::cout<<"[BasisPursuit] Start."<<std::endl;
//		timer.start();
//		kl1p::TBasisPursuitSolver<klab::DoubleReal> bp(tolerance);
//		bp.solve(y, A, x);
//		timer.stop();
//		std::cout<<"[BasisPursuit] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<bp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "BasisPursuit-Signal.csv");	// Write solution to a file.
//
//		// Compute OMP.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[OMP] Start."<<std::endl;
//		timer.start();
//		kl1p::TOMPSolver<klab::DoubleReal> omp(tolerance);
//		omp.solve(y, A, k, x);
//		timer.stop();
//		std::cout<<"[OMP] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<omp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "OMP-Signal.csv");	// Write solution to a file.
//
//		// Compute ROMP.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[ROMP] Start."<<std::endl;
//		timer.start();
//		kl1p::TROMPSolver<klab::DoubleReal> romp(tolerance);
//		romp.solve(y, A, k, x);
//		timer.stop();
//		std::cout<<"[ROMP] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<romp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "ROMP-Signal.csv");	// Write solution to a file.
//
//		// Compute CoSaMP.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[CoSaMP] Start."<<std::endl;
//		timer.start();
//		kl1p::TCoSaMPSolver<klab::DoubleReal> cosamp(tolerance);
//		cosamp.solve(y, A, k, x);
//		timer.stop();
//		std::cout<<"[CoSaMP] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<cosamp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "CoSaMP-Signal.csv");	// Write solution to a file.
//
//		// Compute Subspace-Pursuit.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[SubspacePursuit] Start."<<std::endl;
//		timer.start();
//		kl1p::TSubspacePursuitSolver<klab::DoubleReal> sp(tolerance);
//		sp.solve(y, A, k, x);
//		timer.stop();
//		std::cout<<"[SubspacePursuit] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<sp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "SubspacePursuit-Signal.csv");	// Write solution to a file.
//
//		// Compute SL0.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[SL0] Start."<<std::endl;
//		timer.start();
//		kl1p::TSL0Solver<klab::DoubleReal> sl0(tolerance);
//		sl0.solve(y, A, x);
//		timer.stop();
//		std::cout<<"[SL0] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<sl0.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "SL0-Signal.csv");	// Write solution to a file.
//
//		// Compute AMP.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[AMP] Start."<<std::endl;
//		timer.start();
//		kl1p::TAMPSolver<klab::DoubleReal> amp(tolerance);
//		amp.solve(y, A, x);
//		timer.stop();
//		std::cout<<"[AMP] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<amp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "AMP-Signal.csv");	// Write solution to a file.
//
//		// Compute EMBP.
//		std::cout<<"------------------------------"<<std::endl;
//		std::cout<<"[EMBP] Start."<<std::endl;
//		timer.start();
//		kl1p::TEMBPSolver<klab::DoubleReal> embp(tolerance);
//		embp.enableHomogeneous(true);
//		embp.solve(y, A, k, x);
//		timer.stop();
//		std::cout<<"[EMBP] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
//			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
//				  <<"Iterations="<<embp.iterations()<<std::endl;
//		if(bWrite)
//			kl1p::WriteToCSVFile(x, "EMBP-Signal.csv");	// Write solution to a file.
//
//
//		std::cout<<"------------------------------"<<std::endl;
//
//		std::cout<<std::endl;
//		std::cout<<"End of example."<<std::endl;
//	}
//	catch(klab::KException& e)
//	{
//		std::cout<<"ERROR! KLab exception : "<<klab::FormatExceptionToString(e)<<std::endl;
//	}
//	catch(std::exception& e)
//	{
//		std::cout<<"ERROR! Standard exception : "<<klab::FormatExceptionToString(e)<<std::endl;
	}
	catch(...)
	{
		std::cout<<"ERROR! Unknown exception !"<<std::endl;
	}
}

// ---------------------------------------------------------------------------------------------------- //
