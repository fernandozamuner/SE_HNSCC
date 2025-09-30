require 'optparse'
require_relative 'fisher_table'
Signal.trap('PIPE', 'EXIT')
PRECISION = 2

prefix = ''
OptionParser.new{|opts|
  opts.on('--prefix NAME', 'column names prefix') {|value|
    prefix = "#{value}:"
  }
}.parse!(ARGV)
# indices = ARGV[0].split(',').map(&:to_i).map{|idx| idx - 1 }

header = $stdin.readline.chomp.split("\t")

# ARGV: any_columns	A_pos	B_pos	A_neg	B_neg
col_names = ARGV
indices = col_names.map{|col_name| header.index(col_name) }

puts [*col_names[0...-4], "#{prefix}:numSamples", "#{prefix}:ratio", "#{prefix}:neglog10(corr.significance)"].join("\t")
# puts [*header, 'numSamples', 'ratio', '-log10(significance)'].join("\t")

rows = $stdin.readlines
num_samples = rows.size
rows.map{|l|
  row = l.chomp.split("\t")
  *rest, a_pos, b_pos, a_neg, b_neg = row.values_at(*indices)
  a_pos, b_pos, a_neg, b_neg = [a_pos, b_pos, a_neg, b_neg].map(&:to_i)
  ft = FisherTable.by_two_classes(class_a_positive: a_pos, class_b_positive: b_pos, class_a_negative: a_neg, class_b_negative: b_neg)
  corrected_significance = ft.significance * num_samples
  {row: row, rest_columns: rest, fisher_table: ft, total: ft.total, ratio: ft.a_to_b_positive_rate_ratio, logpvalue: -Math.log10(corrected_significance)}
}.sort_by{|info|
  info[:logpvalue]
}.reverse.each{|info|
  puts [*info[:rest_columns], info[:total], info[:ratio]&.round(PRECISION), info[:logpvalue]&.round(PRECISION)].join("\t")
}
