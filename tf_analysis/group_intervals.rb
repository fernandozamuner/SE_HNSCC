$stdin.each_line.map{|l| l.chomp.split("\t") }.chunk{|r|
  r.first(4)
}.map{|grp, rows|
  intervals = rows.map{|r|
    chr, from, to = r.drop(4)
    [chr, from, to]
  }.reject{|chr, from, to|
    chr == '.' && from == '-1' && to == '-1'
  }.map{|chr, from, to|
    "#{chr}:#{from}-#{to}"
  }.join(",")
  info = [*grp, intervals.empty? ? '.' : intervals]
  puts info.join("\t")
}
