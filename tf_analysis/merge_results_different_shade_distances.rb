filename_1 = ARGV[0]
filename_2 = ARGV[1]
def load_significances(filename)
  File.readlines(filename).drop(1).map{|l|
    tf, num_samples, ratio, neglog_significance = l.chomp.split("\t")
    {
      tf: tf,
      num_samples: Integer(num_samples),
      ratio: ratio.empty? ? nil : Float(ratio),
      neglog_significance: neglog_significance.empty? ? nil : Float(neglog_significance),
    }
  }.sort_by{|info| info[:tf] }
end

tbl_1 = load_significances(filename_1)
tbl_2 = load_significances(filename_2)

tfs_1 = tbl_1.map{|info| info[:tf] }
tfs_2 = tbl_2.map{|info| info[:tf] }

raise  unless tfs_1 == tfs_2

puts ['TF', 'numSamples', 'ratio', '-log10(significance)'].join("\t")
tbl_1.zip(tbl_2).map{|info_1, info_2|
  tf = info_1[:tf]
  num_samples = [info_1[:num_samples], info_2[:num_samples]].min
  neglog_significance = [info_1[:neglog_significance], info_2[:neglog_significance]].min
  
  if info_1[:ratio] && info_2[:ratio]
    if info_1[:ratio] <= 1 && info_2[:ratio] <= 1
      ratio = [info_1[:ratio], info_2[:ratio]].max
    elsif info_1[:ratio] >= 1 && info_2[:ratio] >= 1
      ratio = [info_1[:ratio], info_2[:ratio]].min
    else
    end
  else
    ratio = nil
  end
  
  {tf: tf, num_samples: num_samples, ratio: ratio, neglog_significance: neglog_significance}
}.sort_by{|info|
  -info[:neglog_significance]
}.each{|info|
  row = info.values_at(:tf, :num_samples, :ratio, :neglog_significance)
  puts row.join("\t")
}
