#include <Rcpp.h>

#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
void test() { std::cout << "hi" << std::endl; }