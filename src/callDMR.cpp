#include <Rcpp.h>
using namespace Rcpp;

// Forward declaration - remove default parameters
IntegerMatrix findRegionCpp(CharacterVector chr, NumericVector pos, double sep);
DataFrame findBumpsCpp(CharacterVector chr,
                       NumericVector pos,
                       NumericVector scores,
                       double cutoff,
                       double sep,
                       double dis_merge,
                       double pct_sig,
                       int minCG);

//// [[Rcpp::export]] 
IntegerMatrix findRegionCpp(CharacterVector chr,
                            NumericVector pos,
                            double sep) { // = 1000 (default)
  // Check input parameters
  if (chr.size() == 0 || pos.size() == 0 || chr.size() != pos.size()) {
    Rcpp::stop("Input vectors 'chr' and 'pos' must have the same non-zero length.");
  }

  int n = pos.size();
  std::vector<int> start_indices, end_indices;

  // Initialize the start of the first region
  start_indices.push_back(0);

  // Loop through each position to find region split points
  for (int i = 1; i < n; ++i) {
    // If the difference between the current and previous position is greater than the threshold or chromosome changes, it's a new region
    if (chr[i] != chr[i - 1] || std::abs(pos[i] - pos[i - 1]) > sep) {
      end_indices.push_back(i - 1);    // End of current region
      start_indices.push_back(i);     // Start of new region
    }
  }

  // End of the last region
  end_indices.push_back(n - 1);

  // Convert results to matrix form
  int region_count = start_indices.size();
  IntegerMatrix regions(2, region_count);

  for (int i = 0; i < region_count; ++i) {
    regions(0, i) = start_indices[i]; // Start index
    regions(1, i) = end_indices[i];   // End index
  }

  return regions;
}

// [[Rcpp::export]] don nont export
DataFrame findBumpsCpp(CharacterVector chr,
                       NumericVector pos,
                       NumericVector scores,
                       double cutoff,
                       double sep , // = 1000 (default)
                       double dis_merge, //  = 200
                       double pct_sig, //  = 0.3
                       int minCG) { // = 3

  int n = chr.length();
  if (n == 0) {
    stop("Input vectors are empty.");
  }

  LogicalVector flag(n, false);

  // Set flag values
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(scores[i])) {
      flag[i] = (scores[i] < cutoff);
    }
  }

  // If there are no significant positions, return empty DataFrame
  if (sum(flag) == 0) {
    return DataFrame::create();
  }

  // Find regions
  IntegerMatrix regions = findRegionCpp(chr, pos, sep);
  int ncol = regions.ncol();

  // Pre-allocate result containers
  const int init_capacity = 100000; // Initial capacity
  CharacterVector result_chr(init_capacity);
  NumericVector result_start(init_capacity);
  NumericVector result_end(init_capacity);
  NumericVector result_length(init_capacity);
  IntegerVector result_idx_start(init_capacity);
  IntegerVector result_idx_end(init_capacity);

  int result_idx = 0;

  // Loop through each region
  for (int i = 0; i < ncol; ++i) {
    int region_start = regions(0, i);
    int region_end = regions(1, i);

    // Check if there are enough CGs in the region
    if ((region_end - region_start + 1) <= minCG) { // default < minCG
      continue;
    }

    // Filter region based on flag
    LogicalVector flag_region = flag[seq(region_start, region_end)];
    NumericVector pos_region = pos[seq(region_start, region_end)];
    NumericVector scores_region = scores[seq(region_start, region_end)];
    int nn = flag_region.length();

    // Find start and end positions
    std::vector<int> startidx, endidx;

    // Find start positions
    for (int j = 0; j < nn; ++j) {
      if (flag_region[j] && (j == 0 || !flag_region[j - 1])) {
        startidx.push_back(j);
      }
    }

    // Find end positions
    for (int j = 0; j < nn; ++j) {
      if (flag_region[j] && (j == nn - 1 || !flag_region[j + 1])) {
        endidx.push_back(j);
      }
    }

    // Merge adjacent regions
    std::vector<int> merged_startidx, merged_endidx;
    if (!startidx.empty()) {
      merged_startidx.push_back(startidx[0]);
      for (size_t j = 1; j < startidx.size(); ++j) {
        if (pos_region[startidx[j]] - pos_region[endidx[j - 1]] > dis_merge) {
          merged_endidx.push_back(endidx[j - 1]);
          merged_startidx.push_back(startidx[j]);
        }
      }
      merged_endidx.push_back(endidx.back());
    }

    // Check significant percentage
    for (size_t j = 0; j < merged_startidx.size(); ++j) {
      int start = merged_startidx[j];
      int end = merged_endidx[j];
      int total = end - start + 1;

      // Count significant points
      int count = 0;
      for (int k = start; k <= end; ++k) {
        if (scores_region[k] < cutoff) {
          ++count;
        }
      }

      // Check significant point ratio
      if (static_cast<double>(count) / total > pct_sig) {
        if (result_idx == result_chr.size()) {
          // Dynamic reallocation
          int new_capacity = result_chr.size() * 2;
          result_chr = CharacterVector(new_capacity);
          result_start = NumericVector(new_capacity);
          result_end = NumericVector(new_capacity);
          result_length = NumericVector(new_capacity);
          result_idx_start = IntegerVector(new_capacity);
          result_idx_end = IntegerVector(new_capacity);
        }

        // Save results
        result_chr[result_idx] = chr[region_start + start];
        result_start[result_idx] = pos_region[start];
        result_end[result_idx] = pos_region[end];
        result_length[result_idx] = pos_region[end] - pos_region[start] + 1;
        result_idx_start[result_idx] = region_start + start;
        result_idx_end[result_idx] = region_start + end;
        ++result_idx;
      }
    }
  }

  // If no results, return empty DataFrame
  if (result_idx == 0) {
    return DataFrame::create();
  }

  // Return results, trim extra parts
  return DataFrame::create(
    _["chr"] = result_chr[Range(0, result_idx - 1)],
                         _["start"] = result_start[Range(0, result_idx - 1)],
                                                  _["end"] = result_end[Range(0, result_idx - 1)],
                                                                       _["length"] = result_length[Range(0, result_idx - 1)],
                                                                                                  _["idx.start.global"] = result_idx_start[Range(0, result_idx - 1)],
                                                                                                                                          _["idx.end.global"] = result_idx_end[Range(0, result_idx - 1)]
  );
}

// [[Rcpp::export]]
DataFrame calldmrs_turbo(DataFrame DMLresult, //callDMRCpp
                     double p_threshold = 1e-5,
                     int minlen = 50,
                     int minCG = 3,
                     double dis_merge = 100,
                     double pct_sig = 0.5,
                     double sep = 5000) {

  const double delta = 0;

  // Check and remove NA values
  NumericVector stat = DMLresult["stat"];
  LogicalVector not_na = !is_na(stat);
  double na_proportion = mean(not_na);

  //Rcout << "Proportion of non-NA values: " << na_proportion << std::endl; // Necessary output message
  if (na_proportion < 1.0) {
    // Correctly handle the subset of the DataFrame
    DMLresult = DMLresult[not_na];
  }

  // Check if it is a multifactor design
  bool flag_multifactor = false;
  CharacterVector class_attr = DMLresult.attr("class");
  if (class_attr.length() > 0) {
    for (int i = 0; i < class_attr.length(); i++) {
      if (as<std::string>(class_attr[i]) == "DMLtest.multiFactor") {
        flag_multifactor = true;
        break;
      }
    }
  }

  if (dis_merge > minlen) {
    dis_merge = minlen;
  }

  // Calculate scores
  NumericVector scores;
  if (delta > 0) {
    if (flag_multifactor) {
      stop("The test results is based on multifactor design, 'delta' is not supported");
    }

    NumericVector diff = DMLresult["stat"];
    NumericVector diff_se = DMLresult["diff.se"];
    int n = diff.length();
    NumericVector postprob(n);
    scores = NumericVector(n);

    // Vectorized operation for calculating scores
    for (int i = 0; i < n; i++) {
      if (R_IsNA(diff[i]) || R_IsNA(diff_se[i])) {
        scores[i] = NA_REAL;
        continue;
      }
      double p1 = R::pnorm(diff[i] - delta, 0, diff_se[i], 1, 0);
      double p2 = R::pnorm(diff[i] + delta, 0, diff_se[i], 0, 0);
      postprob[i] = p1 + p2;
      scores[i] = 1 - postprob[i];
    }
  } else {
    scores = DMLresult["pval"];
  }

  // Find DMRs - using findBumpsCpp function to process input data
  DataFrame dmrs = findBumpsCpp(DMLresult["chr"],
                                DMLresult["pos"],
                                         scores,
                                         p_threshold,
                                         sep,           // 5000 Fixed value for sep parameter
                                         dis_merge,      // Merge distance
                                         pct_sig,        // Significance percentage
                                         minCG);         // Minimum number of CpGs



  // Print processing progress and results
  // Rcout << "Calling DMRs  ..." << std::endl;
  // Rcout << "Number of DMRs found: " << dmrs.nrows() << std::endl;

  // Check if any DMRs were found
  if (dmrs.nrows() == 0) {
    Rcout << "No DMRs found." << std::endl;
    return DataFrame::create();  // Return an empty data frame
  }

  // Extract necessary columns from dmrs
  NumericVector lengths = dmrs["length"];        // Length of DMRs
  IntegerVector idx_start = dmrs["idx.start.global"];  // Starting position index
  IntegerVector idx_end = dmrs["idx.end.global"];      // Ending position index

  // Check if the lengths of start and end indices match
  if (idx_start.length() != idx_end.length()) {
    stop("Error: idx_start and idx_end must be of the same length.");
  }

  // Calculate the number of CpG sites in each region
  IntegerVector nCG = idx_end - idx_start + 1;

  // Check the size of lengths and nCG vectors to ensure they're valid
  if (lengths.length() != nCG.length()) {
    stop("Error: lengths and nCG vectors must have the same length.");
  }

  // Filter the rows based on lengths and CpG count
  LogicalVector ix_good = (lengths > minlen) & (nCG >= minCG); // default nCG > minCG

  // Create new vectors with filtered data
  CharacterVector filtered_chr = as<CharacterVector>(dmrs["chr"]);
  NumericVector filtered_start = as<NumericVector>(dmrs["start"]);
  NumericVector filtered_end = as<NumericVector>(dmrs["end"]);
  NumericVector filtered_lengths = lengths;
  IntegerVector filtered_idx_start = idx_start;
  IntegerVector filtered_idx_end = idx_end;

  // Apply the filter
  filtered_chr = filtered_chr[ix_good];
  filtered_start = filtered_start[ix_good];
  filtered_end = filtered_end[ix_good];
  filtered_lengths = filtered_lengths[ix_good];
  filtered_idx_start = filtered_idx_start[ix_good];
  filtered_idx_end = filtered_idx_end[ix_good];

  // Create the filtered DataFrame
  dmrs = DataFrame::create(
    _["chr"] = filtered_chr,
    _["start"] = filtered_start,
    _["end"] = filtered_end,
    _["length"] = filtered_lengths,
    _["idx.start.global"] = filtered_idx_start,
    _["idx.end.global"] = filtered_idx_end
  );

  // Print information on the filtering process
  // Rcout << "After filtering - Number of rows: " << dmrs.nrows() << std::endl;
  //Rcout << "Number of DMRs found: " << dmrs.nrows() << std::endl;

  // Extract necessary columns from dmrs
  IntegerVector filtered_nCG = filtered_idx_end - filtered_idx_start + 1;

  // Print filtered information
  // Rcout << "Filtered to " << dmrs.nrows() << " rows" << std::endl;
  // Rcout << "\nDetailed index information for each DMR:" << std::endl;

  // for(int i = 0; i < filtered_idx_start.length(); i++) {
  //   Rcout << "DMR " << i + 1 << ": ";
  //   Rcout << "Start=" << filtered_idx_start[i] << ", ";
  //   Rcout << "End=" << filtered_idx_end[i] << ", ";
  //   Rcout << "CpG count=" << filtered_nCG[i] << std::endl;
  // }


  // 6. Subsequent code uses dmrs instead of dmrs
  if (flag_multifactor) {
    NumericVector areaStat(dmrs.nrows());
    NumericVector stat = DMLresult["stat"];

    for (int i = 0; i < dmrs.nrows(); i++) {
      double sum_stat = 0;
      for (int j = filtered_idx_start[i]; j <= filtered_idx_end[i]; j++) {
        if (!R_IsNA(stat[j])) {
          sum_stat += stat[j];
        }
      }
      areaStat[i] = sum_stat;
    }

    // Return results
    return DataFrame::create(
      _["chr"] = dmrs["chr"],
                     _["start"] = dmrs["start"],
                                      _["end"] = dmrs["end"],
                                                     _["length"] = dmrs["length"],
                                                                       _["nCG"] = filtered_nCG,  // Use newly computed nCG
                                                                       _["areaStat"] = areaStat
    );
  } else {
    NumericVector meanMethy1(dmrs.nrows());
    NumericVector meanMethy2(dmrs.nrows());
    NumericVector areaStat(dmrs.nrows());
    NumericVector mu1 = DMLresult["mu1"];
    NumericVector mu2 = DMLresult["mu2"];
    NumericVector stat = DMLresult["stat"];

    for (int i = 0; i < dmrs.nrows(); i++) {
      double sum1 = 0, sum2 = 0, sum_stat = 0;
      int count = 0;

      // Print current DMR's range
      // Rcout << "\nProcessing DMR " << i + 1 << "/" << dmrs.nrows() << ":" << std::endl;
      // Rcout << "Start index: " << filtered_idx_start[i] << ", End index: " << filtered_idx_end[i] << std::endl;

      for (int j = filtered_idx_start[i]; j <= filtered_idx_end[i]; j++) {
        if (j >= mu1.length() || j >= mu2.length() || j >= stat.length()) {
          stop("Index out of bounds error in methylation calculation");
        }
        // Print the current position's value
        // Rcout << "Position " << j << ": ";
        // Rcout << "mu1=" << mu1[j] << ", mu2=" << mu2[j] << ", stat=" << stat[j];

        if (!R_IsNA(mu1[j]) && !R_IsNA(mu2[j])) {
          sum1 += mu1[j];
          sum2 += mu2[j];
          sum_stat += stat[j];
          count++;
          //Rcout << " (included)";
        //} else {
          //Rcout << " (NA - skipped)";
        }
        //Rcout << std::endl;
      }

      // Calculate mean methylation and areaStat
      meanMethy1[i] = (count > 0) ? sum1 / count : NA_REAL;
      meanMethy2[i] = (count > 0) ? sum2 / count : NA_REAL;
      areaStat[i] = sum_stat;

      // Rcout << "DMR " << i + 1 << " summary:" << std::endl;
      // Rcout << "  Valid CpG sites: " << count << std::endl;
      // Rcout << "  Mean methylation 1: " << meanMethy1[i] << std::endl;
      // Rcout << "  Mean methylation 2: " << meanMethy2[i] << std::endl;
      //Rcout << "  Area statistic: " << areaStat[i] << std::endl;
      //Rcout << "----------------------------------------" << std::endl;
    }

    // Return results
    return DataFrame::create(
      _["chr"] = dmrs["chr"],
                     _["start"] = dmrs["start"],
                                      _["end"] = dmrs["end"],
                                                     _["length"] = dmrs["length"],
                                                                       _["nCG"] = filtered_nCG,  // Use newly computed nCG
                                                                       _["meanMethy1"] = meanMethy1,
                                                                       _["meanMethy2"] = meanMethy2,
                                                                       _["diff.Methy"] = meanMethy1 - meanMethy2,
                                                                       _["areaStat"] = areaStat
    );
  }
}
