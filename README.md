# Nonlinear-Optimization-Algorithms
MATLAB implementations of various nonlinear programming algorithms.

---

This repository contains MATLAB implementations of a variety of popular nonlinear programming algorithms, many of which can be found in *Numerical Optimization* by Nocedal and Wright, a text that I highly recommend.

List of algorithms implemented:

1. [line-search](https://en.wikipedia.org/wiki/Line_search) ([simple Wolfe, strong Wolfe](https://en.wikipedia.org/wiki/Wolfe_conditions), [Mor&#233;-Thuente](http://dl.acm.org/citation.cfm?id=192132))
2. [steepest descent](https://en.wikipedia.org/wiki/Method_of_steepest_descent)
3. [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization)
4. [Dogleg method](http://www.numerical.rl.ac.uk/people/nimg/course/lectures/raphael/lectures/lec7slides.pdf)
5. [Steihaug-Toint conjugate gradient trust region method](http://www.numerical.rl.ac.uk/people/nimg/course/lectures/raphael/lectures/lec7slides.pdf)
6. [BFGS](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm)
7. [limited-memory BFGS](https://en.wikipedia.org/wiki/Limited-memory_BFGS)
8. [Gauss-Newton method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm)

All of the algorithms are heavily commented (possibly to a fault), but I wanted someone in the midst of a nonlinear programming class to be able to read through the code and understand it decently well.  Although I have done my best to implement these algorithms with efficiency in mind (within the confines of MATLAB's inherent deficiencies in this regard), this repository is far more valuable as a teaching tool than as a performance-centric library.

Due to the algorithms being so heavily commented, many implementation details are contained within the code as comments instead of in a README.

Some day, I will include a demo folder that demonstrates the correctness and performance of each algorithm on a set of representative problems, and I will create a README with implementation details for each algorithm, to be located in the src folder.

Some day! :)
