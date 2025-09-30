require 'shellwords'
Dir.glob("source_data/hg19_cismotifs/*.bed").group_by{|fn|
	File.basename(fn).split(".").first
}.map{|tf, fns|
  abc_fns = fns.select{|fn|
    ['A', 'B', 'C'].include?(File.basename(fn).split(".")[1])
  }
  [tf, abc_fns]
}.reject{|tf,fns|
  fns.empty?
}.each{|tf, fns|
  cmd = [
    "cat #{fns.shelljoin}",
    'bedtools sort',
    'ruby filter_bed_by_length.rb 0 10000',
    'bedtools merge',
    'ruby filter_bed_by_length.rb 50 10000'
  ].join(' | ')
  output_fn = "results/cistrome_hg19_abc/#{tf}.bed"
  puts "#{cmd} > #{output_fn}"
}
