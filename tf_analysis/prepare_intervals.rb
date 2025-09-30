require 'fileutils'
require_relative 'shade_generator'

def load_regions(filename)
  File.open(filename){|f|
    header = f.readline.chomp.split("\t")
    columns = ['seqnames', 'start', 'end', 'sample_type', 'name', 'ID']
    fields =  [:chr, :from, :to, :sample_type, :region_type, :id]
    columns_indices = columns.map{|col_name| header.index(col_name) }
    f.each_line.map{|l|
      row = l.chomp.split("\t")
      res = fields.zip( row.values_at(*columns_indices) ).to_h
      res[:from] = Integer(res[:from])
      res[:to] = Integer(res[:to])
      res
    }
  }
end

###########

intervals_fns = ARGV # ['source_data/samples_3456_lilly_results_with_DE.annotated.tsv', 'source_data/cell_lines_lilly_results.tsv']
chrom_sizes = load_chrom_sizes('source_data/hg19.chrom.sizes')
all_regions = intervals_fns.flat_map{|fn| load_regions(fn) }
FileUtils.mkdir_p 'intervals'

SAMPLE_TYPES = {
  'CellLine': ['CellLine047', 'CellLine090'],
  'CellLine047': ['CellLine047'],
  'CellLine090': ['CellLine090'],
  'tumor': ['tumor'],
  'UPPP': ['UPPP']
}

["SE", "enhancer", "promoter"].each{|region_type|
  SAMPLE_TYPES.each{|resulting_sample_type, sample_type_options|
    regions = all_regions.select{|info|
      sample_type_options.include?(info[:sample_type]) && (info[:region_type] == region_type)
    }
    kwargs = {sample_type: resulting_sample_type, region_type: region_type}
    store_regions_and_shades(regions, "intervals/#{region_type}_#{resulting_sample_type}", chrom_sizes: chrom_sizes, **kwargs)
  }
}
