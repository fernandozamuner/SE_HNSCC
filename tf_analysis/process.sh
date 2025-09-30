# Merge cistrome ABC-categories for each factor and drop too long (and too short 50-10000 nt are allowed) peaks
ruby cistrome_merger_cmd.rb | bash

ruby prepare_intervals.rb source_data/samples_3456_lilly_results_with_DE.annotated.tsv source_data/cell_lines_lilly_results.tsv
ruby prepare_dse_intervals.rb

ruby count_peaks_covered.rb intervals/*.bed > results/intervals_overlap_counts.tsv
ruby all_significances.rb  results/intervals_overlap_counts.tsv

ruby count_peaks_covered.rb dse_intervals/*.bed > results/dse_intervals_overlap_counts.tsv

ruby calc_significance_sites_vs_nosites.rb  results/intervals_overlap_counts.tsv  results
ruby calc_significance_all_vs_all.rb  results/intervals_overlap_counts.tsv  results

ruby calc_significance_sites_vs_nosites.rb  results/dse_intervals_overlap_counts.tsv  results_dse

ruby prepare_intervals_new.rb intervals_new
mkdir -p results_new
ruby count_peaks_covered.rb intervals_new/*.bed > results_new/intervals_overlap_counts.tsv
ruby calc_significance_sites_vs_nosites.rb  results_new/intervals_overlap_counts.tsv  results_new


ruby prepare_intervals_new.rb intervals_bySource
mkdir -p results_bySource
ruby count_peaks_covered.rb intervals_bySource/*.bed > results_bySource/intervals_overlap_counts.tsv
ruby calc_significance_sites_vs_nosites.rb  results_bySource/intervals_overlap_counts.tsv  results_bySource
