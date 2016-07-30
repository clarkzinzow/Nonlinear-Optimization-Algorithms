# Nonlinear-Optimization-Algorithms
MATLAB implementations of a variety of nonlinear programming algorithms.

---

This repository contains MATLAB implementations of a variety of popular nonlinear programming algorithms, many of which can be found in *Numerical Optimization" by Nocedal and Wright, a text that I highly recommend.

List of algorithms implemented:

1. line-search (simple Wolfe, strong Wolfe, Mor&#233;-Thuente)
2. steepest descent
3. Newton's method
4. Dogleg method
5. Steihaug-Toint conjugate gradient trust region method
6. BFGS
7. limited-memory BFGS
8. Gauss-Newton method

All of the algorithms are heavily commented (possibly to a fault), but I wanted someone in the midst of a nonlinear programming class to be able to read through the code and understand it decently well.  Although I have done my best to implement these algorithms with efficiency in mind (within the confines of MATLAB's inherent deficiencies in this regard), this repository are far more valuable as a teaching tool than a performance-centric library.

In the near future, I will include a demo folder that demonstrates the correctness and performance of each algorithm on a set of representative problems.