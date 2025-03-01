// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calldmrs_turbo
DataFrame calldmrs_turbo(DataFrame DMLresult, double p_threshold, int minlen, int minCG, double dis_merge, double pct_sig, double sep);
RcppExport SEXP _methylTracer_calldmrs_turbo(SEXP DMLresultSEXP, SEXP p_thresholdSEXP, SEXP minlenSEXP, SEXP minCGSEXP, SEXP dis_mergeSEXP, SEXP pct_sigSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type DMLresult(DMLresultSEXP);
    Rcpp::traits::input_parameter< double >::type p_threshold(p_thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type minlen(minlenSEXP);
    Rcpp::traits::input_parameter< int >::type minCG(minCGSEXP);
    Rcpp::traits::input_parameter< double >::type dis_merge(dis_mergeSEXP);
    Rcpp::traits::input_parameter< double >::type pct_sig(pct_sigSEXP);
    Rcpp::traits::input_parameter< double >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(calldmrs_turbo(DMLresult, p_threshold, minlen, minCG, dis_merge, pct_sig, sep));
    return rcpp_result_gen;
END_RCPP
}
// computeStatCpp
DataFrame computeStatCpp(NumericVector mean_1, NumericVector mean_2, NumericVector var1, NumericVector var2);
RcppExport SEXP _methylTracer_computeStatCpp(SEXP mean_1SEXP, SEXP mean_2SEXP, SEXP var1SEXP, SEXP var2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mean_1(mean_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean_2(mean_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var1(var1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var2(var2SEXP);
    rcpp_result_gen = Rcpp::wrap(computeStatCpp(mean_1, mean_2, var1, var2));
    return rcpp_result_gen;
END_RCPP
}
// fill_missing_with_mean
IntegerMatrix fill_missing_with_mean(IntegerMatrix chunk, NumericVector chunk_mean);
RcppExport SEXP _methylTracer_fill_missing_with_mean(SEXP chunkSEXP, SEXP chunk_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type chunk(chunkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chunk_mean(chunk_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(fill_missing_with_mean(chunk, chunk_mean));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_methylTracer_calldmrs_turbo", (DL_FUNC) &_methylTracer_calldmrs_turbo, 7},
    {"_methylTracer_computeStatCpp", (DL_FUNC) &_methylTracer_computeStatCpp, 4},
    {"_methylTracer_fill_missing_with_mean", (DL_FUNC) &_methylTracer_fill_missing_with_mean, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_methylTracer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
