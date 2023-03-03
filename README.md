# Bayesian Sketches 

This project provides a codebase with which the Count-[1], CountMin-[2], CB-, CCA-, CCB-, Seq-[3], and PR-Sketches[4] can be compared to each other in a simulation environment.
The CB-, CCA-, and CCB-Sketches are the algorithms proposed by the paper[9] associated to this project.

## What is inside

This project can be used to read in real traces (until now, only CAIDA-formatted traces are supported, that is, the datastream is a sequence of rows, each of which is formatted as [timestamp, ID, payload]), or create synthetic traces. 

After having defined the trace to be used, the trace can be fed to a compressing-datastructure, of which there are three kinds implemented: The simple list of counter-arrays as used by CountMin-, PR-, CB-, CCA-, and CCB-Sketches. And the randomly-normalizing modified list of counter arrays as used by the Seq-Sketch. The two-level-hashing modified list of counter arrays as used by the Count-Sketch. All compressive datastructures, except for the latter used by the Count-Sketch, also keep track of the number of distinct keys that have been associated to each counter, or in other words it keeps track of cardinality information. One can adjust the accuracy of the cardinality information by adjusting the noise-parameter which specifies the fraction of distinct keys that do not get observed by the procedure and hence yield faulty cardinality information. As of right now, this perturbation of cardinality information is not supported for the PR- and Seq-Sketch and does not play a role for the Count- and CountMin-Sketches. 
It is also possible to simulate cardinality counter noise stemming from a bloom filter implementation keeping track of existing IDs. 

After having inserted the trace into the compressing datastructure, one chooses which query-algorithm he or she wishes to use, under the constraint that the compressive datastructure and the query algorithm must be compatible.

Next, an execution engine for the pair of (datastructure, query-algorithm) is created which computes the predicitions of aggregate volume sizes given by the query algorithm. The predictions are then compared with the ground truth and error statistics are computed.

Lastly, error statsitics and metadata can be dumped into a csv file.

The main.cpp file provides and example of how this project is intended to be used.

## What is not inside

This project does not provide highly optimized code for sketching algorithms as it trades off performance for modularity and ease of understanding. The authors attempted to make the playing-field as even as possible by implementing the different query algorithms in a common style and using a common codebase. Seq- and PR-Sketch are two exceptions as they require specialized algorithms to be viable at all and have hence been implemented using a custom compressive sensing library[7] and armadillo[5, 6] respecitvely. Regarding the Seq- and PR-Sketch, it is unclear whether the performance of this implementation matches what has been implemented as part of the respecive publications since time and memory measurements in those publications where either ambigous or missing and furthermore the code attached to the respecive publications was not available at the time of writing. Specifically regarding the PR-Sketch, it is unclear which algorithm should be used to solve the problem especially since it deals with sparse matrices, which reduces the amount of publicly available code fit to the problem. In the PR-Sketch publication it is mentioned that the Eigen library was used to solve the problem, but this does not explain what algorithm has been used.

This project does not provide highly parallelized code. While parallelization is present to a certain degree, it is suboptimal for long traces and will cause a higher memory consumption than expected. For the medium-sized experiments reported in the paper associated to this project, this is not much of an issue.

## Requirements to run the project

The project must first of all be compiled. To that end, besides the basic requirements to compile C++ code, one requires the following external libraries: Eigen 3.4.0 (older versions should also be fine), Armadillo 10.8.2 (some older versions should also be fine), and KL1p 0.4.2. Installing Eigen and Armadillo is fairly straightforward while KL1p is a bit more complicated and requires the user to manually compile the library on their computer. All instructions are given on the websites of the respective libraries.

The C++ dialect is GNU++17 and the compiler used by the authors is Apple clang version 13.0.0 (clang-1300.0.29.30) with arm64-apple-darwin21.2.0 as target architecture. Other compilers and targets should work just as fine if not better.

The project makes heavy use of hashing and hence a MurmurHash3 implementation was used. There is no need to install or download MurmurHash3 as the relevant code is already present in the project as it has been copied from its original online repository[8].

## License

The code of this project, with exception of the MurmurHash3 files, is distributed under the BSD 3-Clause License, which would be:

BSD 3-Clause License

Copyright (c) [2022], [Francesco Da Dalt]

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



## Authors

Francesco Da Dalt is the main author of this project. Of course the MurmurHash3 implementation and the other libraries used are external sources and have been written by different authors as is shown by the references in this document.


# References

[1] Charikar, M., Chen, K., Farach-Colton, M. (2002). Finding Frequent Items in Data Streams

[2] Graham Cormode and S. Muthukrishnan. (2005). An improved data stream summary: the count-min sketch and its applications

[3] Huang, Q., Sheng, S., Chen, X., Bao, Y., Zhang, R., Xu, Y., and Zhang, G. (2021). Toward Nearly-Zero-Error Sketching via Compressive Sensing

[4] Siyuan Sheng, Qun Huang, Sa Wang, and Yungang Bao. (2021). PR-sketch: monitoring per-key aggregation of streaming data with nearly full accuracy

[5] Conrad Sanderson and Ryan Curtin. (2018). A User-Friendly Hybrid Sparse Matrix Class in C++

[6] Conrad Sanderson and Ryan Curtin. (2016). Armadillo: a template-based C++ library for linear algebra

[7] René Gebel. (2012). KL1p: A portable C++ library for Compressed Sensing

[8] Austin Appleby. (2016). https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp

[9] Francesco Da Dalt, Simon Scherrer, Adrian Perrig. (2022). Bayesian Sketches for Volume Estimation in Data Streams. https://doi.org/10.14778/3574245.3574252

