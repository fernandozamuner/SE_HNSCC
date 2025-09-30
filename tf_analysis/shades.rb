header = ARGV.delete('--header')

gap_size = Integer(ARGV[0])
chrom_sizes = nil
if ARGV.size == 2
  genome_sizes_fn = ARGV[1]
  chrom_sizes = File.readlines(genome_sizes_fn).map{|l|
    chr, sz = l.chomp.split("\t")
    [chr, sz.to_i]
  }.to_h
end

puts $stdin.readline  if header

$stdin.each_line{|l|
  row = l.chomp.split("\t")
  chr = row[0]
  from = row[1].to_i
  to = row[2].to_i
  rest = row.drop(3)
  interval_length = to - from

  left_flank = [from - gap_size - interval_length, from - gap_size]
  right_flank = [to + gap_size, to + gap_size + interval_length]
  if chrom_sizes
    left_flank = left_flank.map{|coord| [coord, 0].max }
    right_flank = right_flank.map{|coord| [coord, chrom_sizes[chr]].min }
  end

  puts [chr, *left_flank, *rest].join("\t")
  puts [chr, *right_flank, *rest].join("\t")
}
