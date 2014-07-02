#include <Rcpp.h>
using namespace Rcpp;

// Given a sorted pvalue, z and lfdr vectors, smooth to force the
// lfdr to be increasing with increasing p-value
//
// In any case where pvalue[j] < pvalue[i], abs(z[i] - z[j]) < window,
// and lfdr[j] > lfdr[j], replace lfdr[i] with lfdr[j].
//
// Note that the input *must* be sorted so that pvalue is increasing
// [[Rcpp::export]]
NumericVector monoSmooth(NumericVector pvalue, NumericVector z,
                        NumericVector lfdr, double window) {
    int n = pvalue.size();

    NumericVector out(n);
    
    try {
        // check sizes of vectors match
        if (n != z.size()) {
            throw std::range_error("Lengths of pvalue and z vectors must match");
        }
        if (n != lfdr.size()) {
            throw std::range_error("Lengths of pvalue and lfdr vectors must match");
        }

        double innermax, lbound, ubound;

        for (int i = 0; i < n; i++) {
            innermax = lfdr[i];
            lbound = z[i] - window;
            ubound = z[i] + window;
            for (int j = 0; j < i; j++) {
                if (z[j] > lbound && z[j] < ubound && lfdr[j] > innermax) {
                    innermax = lfdr[j];
                }
            }
            out[i] = innermax;
        }
    } catch (std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }

    return out;
}
