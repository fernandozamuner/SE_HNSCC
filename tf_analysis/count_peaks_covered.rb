enhancers_globs = ARGV
raise 'Specify globs with covering regions (enhancers or control)'  if enhancers_globs.empty?
folder_tf = 'results/cistrome_hg19_abc'
tfs = Dir.glob("#{folder_tf}/*.bed").map{|fn| File.basename(fn, '.bed') }.sort
enhancer_fns = enhancers_globs.flat_map{|glob|
  Dir.glob(glob)
}.sort_by{|fn|
  File.basename(fn, '.bed')
}

enhancer_names = enhancer_fns.map{|fn|
  File.basename(fn, '.bed')
}

header = ['TF', *enhancer_names, *enhancer_names.map{|name| "neg:#{name}"}]
puts header.join("\t")

tfs.each{|tf|
  counts_all_enhancers = enhancer_fns.map{|enhancers_fn|
    cmd = "bedtools intersect -u -a #{enhancers_fn} -b #{folder_tf}/#{tf}.bed | wc -l"
    Integer(`#{cmd}`)
  }
  negcounts_all_enhancers = enhancer_fns.map{|enhancers_fn|
    cmd = "bedtools intersect -v -a #{enhancers_fn} -b #{folder_tf}/#{tf}.bed | wc -l"
    Integer(`#{cmd}`)
  }
  infos = [tf, *counts_all_enhancers, *negcounts_all_enhancers]
  puts infos.join("\t")
  $stderr.print '.'
}
