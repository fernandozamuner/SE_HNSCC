min_length, max_length = ARGV.pop(2).map{|x| Integer(x) }

ARGF.each_line{|l|
  l.chomp!
  chr, s, f, *rest = l.split("\t", 4)
  len = Integer(f) - Integer(s)
  puts l  if (min_length <= len) && (len <= max_length)
}
