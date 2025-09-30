require 'fileutils'
enhancers_globs = ARGV
raise 'Specify globs with covering regions (enhancers or control)'  if enhancers_globs.empty?
folder_tf = 'results/cistrome_hg19_abc'
tfs = Dir.glob("#{folder_tf}/*.bed").map{|fn| File.basename(fn, '.bed') }.sort
enhancer_fns = enhancers_globs.flat_map{|glob|
  Dir.glob(glob)
}.sort_by{|fn|
  File.basename(fn, '.bed')
}

FileUtils.mkdir_p('peaks_covered')
enhancer_fns.each{|enhancers_fn|
  puts enhancers_fn
  enhancers_type = File.basename(enhancers_fn, '.bed')
  common = nil
  column_by_tf = {}
  tfs.each{|tf|
    print('.')
    cmd = "cut -d '\t' -f1-4 #{enhancers_fn} | bedtools intersect -loj -a stdin -b #{folder_tf}/#{tf}.bed | ruby group_intervals.rb" # > peaks_covered/#{enhancers_type}@#{tf}
    tf_intersections = `#{cmd}`.lines.map{|l| l.chomp.split("\t") }
    common ||= tf_intersections.map{|l| l.first(4) }
    column_by_tf[tf] = tf_intersections.map{|l| l.drop(4) }
  }
  File.open("peaks_covered/#{enhancers_type}.bed", 'w'){|fw|
    header = ['chr', 'start', 'stop', 'name', *tfs]
    fw.puts header.join("\t")
    common.zip(*tfs.map{|tf| column_by_tf[tf] }).each{|row|
      fw.puts row.join("\t")
    }
  }
}
