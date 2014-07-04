#include <Rcpp.h>
using namespace Rcpp;

// Given a sorted pvalue, z and lfdr vectors, smooth to force the
// density to be decreasing with increasing p-value
//
// In any case where pvalue[j] < pvalue[i], abs(z[i] - z[j]) < window,
// and density[j] < density[j], replace lfdr[i] with lfdr[j].
//
// Note that the input *must* be sorted so that pvalue is increasing
// [[Rcpp::export]]
NumericVector monoSmooth(NumericVector pvalue, NumericVector z,
                         NumericVector density, double window) {
    int n = pvalue.size();

    try {
        // check sizes of vectors match
        if (n != z.size()) {
            throw std::range_error("Lengths of pvalue and z vectors must match");
        }
        if (n != density.size()) {
            throw std::range_error("Lengths of pvalue and density vectors must match");
        }

        double innermin, lbound, ubound;

        for (int i = 0; i < n; i++) {
            innermin = density[i];
            lbound = z[i] - window;
            ubound = z[i] + window;
            for (int j = 0; j < i; j++) {
                if (z[j] > lbound && z[j] < ubound && density[j] < innermin) {
                    innermin = density[j];
                }
            }
            density[i] = innermin;
        }
    } catch (std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }

    return density;
}
