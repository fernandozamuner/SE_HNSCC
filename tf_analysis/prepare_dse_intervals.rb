require 'fileutils'
require_relative 'shade_generator'

DSEInterval = Struct.new(:log_fc, :log_cpm, :f, :pvalue, :fdr, :id, :chr, :start, :stop, :width, :strand) do
  def self.from_string(str)
    log_fc, log_cpm, f, pvalue, fdr, id, chr, start, stop, width, strand = str.chomp.split(',')
    self.new(Float(log_fc), Float(log_cpm), Float(f), Float(pvalue), Float(fdr), id, chr, Integer(start), Integer(stop), Integer(width), strand)
  end

  def self.all_from_file(fn)
    File.readlines(fn).drop(1).map{|l| self.from_string(l) }
  end

  def as_dict
    {chr: chr, from: start, to: stop, id: id}
  end
end

chrom_sizes = load_chrom_sizes('source_data/hg19.chrom.sizes')

FileUtils.mkdir_p 'dse_intervals'
['../results/DSE_for_TF/P_by_ChipSeq.csv', '../results/DSE_for_TF/E_by_ChipSeq.csv', '../results/DSE_for_TF/SE_by_ChipSeq.csv'].each{|fn|
  regions = DSEInterval.all_from_file(fn).select{|interval| interval.fdr <= 0.05 }
  regions_up = regions.select{|interval| interval.log_fc > 0 }.map(&:as_dict)
  regions_down = regions.select{|interval| interval.log_fc < 0 }.map(&:as_dict)
  basename = File.join("dse_intervals/", File.basename(fn, '.csv'))
  store_regions_and_shades(regions_up, "#{basename}_TUMOR", chrom_sizes: chrom_sizes)
  store_regions_and_shades(regions_down, "#{basename}_UPPP", chrom_sizes: chrom_sizes)
}
