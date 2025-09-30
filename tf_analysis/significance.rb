require_relative 'fisher_table'

MAX_LOGP = 1000 # replace Infinity as a log-significance with a large value (not to screw excel)

def take_only(values)
  raise 'All values should be the same'  unless values.uniq.size == 1
  values.first
end

def merge_info(info_1, info_2)
  num_samples = [info_1[:num_samples], info_2[:num_samples]].min
  neglog_significance = [info_1[:neglog_significance], info_2[:neglog_significance]].compact.min
  a_pos = take_only([info_1, info_2].map{|info| info[:a_pos] })
  a_neg = take_only([info_1, info_2].map{|info| info[:a_neg] })
  
  if info_1[:ratio] && info_2[:ratio]
    if info_1[:ratio] <= 1 && info_2[:ratio] <= 1
      ratio = [info_1[:ratio], info_2[:ratio]].max
    elsif info_1[:ratio] >= 1 && info_2[:ratio] >= 1
      ratio = [info_1[:ratio], info_2[:ratio]].min
    else
      neglog_significance = nil
      ratio = nil
    end
  else
    ratio = nil
  end
  
  {num_samples: num_samples, ratio: ratio, neglog_significance: neglog_significance, a_pos: a_pos, a_neg: a_neg}
end

# column_names: [A+, B+, A-, B-]
def calc_single_tf_significance(row, header, column_names:)
  indices = column_names.map{|col_name| header.index(col_name) }
  a_pos, b_pos, a_neg, b_neg = row.values_at(*indices).map(&:to_i)
  fisher_table = FisherTable.by_two_classes(class_a_positive: a_pos, class_b_positive: b_pos, class_a_negative: a_neg, class_b_negative: b_neg)
  {
    neglog_significance: [-Math.log10(fisher_table.significance), 0].max,
    ratio: fisher_table.a_to_b_odds_ratio,
    num_samples: fisher_table.total,
    a_pos: a_pos, b_pos: b_pos,
    a_neg: a_neg, b_neg: b_neg,
  }
end

def correct_neglog_significance(neglog_significance, bonferroni_correction)
  [[(neglog_significance || 0) - Math.log10(bonferroni_correction), 0].max, MAX_LOGP].min
end
