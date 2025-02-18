#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame computeStatCpp(NumericVector mean_1, 
                        NumericVector mean_2,
                        NumericVector var1,
                        NumericVector var2) {
    
    int n = mean_1.length();
    NumericVector dif(n);
    NumericVector vv(n);
    NumericVector se(n);
    NumericVector stat(n);
    NumericVector pval(n);
    
    // diff
    for(int i = 0; i < n; i++) {
        dif[i] = mean_1[i] - mean_2[i];
        vv[i] = var1[i] + var2[i];
        if(vv[i] < 1e-5) vv[i] = 1e-4;
        if(vv[i] > 1 - 1e-5) vv[i] = 1 - 1e-4;
        se[i] = sqrt(vv[i]);
        stat[i] = dif[i] / se[i];
        pval[i] = 2 * R::pnorm(-std::abs(stat[i]), 0.0, 1.0, 1, 0);
    }
    Environment stats("package:stats");
    Function p_adjust = stats["p.adjust"];
    NumericVector fdr = p_adjust(pval, Named("method") = "fdr");
    return DataFrame::create(
        Named("mu1") = mean_1,
        Named("mu2") = mean_2,
        Named("diff") = dif,
        Named("diff.se") = se,
        Named("stat") = stat,
        Named("pval") = pval,
        Named("fdr") = fdr
    );
}

// [[Rcpp::export]]
IntegerMatrix fill_missing_with_mean(IntegerMatrix chunk, NumericVector chunk_mean) {
  int nrows = chunk.nrow();
  int ncols = chunk.ncol();
  
  for (int i = 0; i < nrows; i++) {
    double mean_value = chunk_mean[i];
    
    for (int j = 0; j < ncols; j++) {
      if (IntegerVector::is_na(chunk(i, j))) {
        chunk(i, j) = static_cast<int>(mean_value);
      }
    }
  }
  
  return chunk;
}
