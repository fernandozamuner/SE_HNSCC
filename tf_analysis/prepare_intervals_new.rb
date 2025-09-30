require 'csv'
require 'fileutils'
require_relative 'shade_generator'

raise "Specify destination folder"  unless destination_folder = ARGV[0]

def load_regions(filename)
  CSV.readlines(filename, headers: true).map(&:to_h).map{|region|
    {
      id: region['ID'],
      chr: region['seqnames'],
      from: Integer(region['start']),
      to: Integer(region['end']),
      logFC: Float(region['logFC']),
      logCPM: Float(region['logCPM']),
      F: Float(region['F']),
      pvalue: Float(region['PValue']),
      fdr: Float(region['FDR']),
      width: Integer(region['width']),
      source: region['source'],
      samples: {
        coverage: {
          UPPP: [
            Integer(region['X3_h3k27ac_bowtie2_hg19_sort.bam']),
            Integer(region['X4_h3k27ac_bowtie2_hg19_sort.bam']),
          ],
          tumor: [
            Integer(region['X5_h3k27ac_bowtie2_hg19_sort.bam']),
            Integer(region['X6_h3k27ac_bowtie2_hg19_sort.bam']),
          ],
        },
        cpm: {
          UPPP: [
            Float(region['X3_h3k27ac_bowtie2_hg19_sort.bam.cpm']),
            Float(region['X4_h3k27ac_bowtie2_hg19_sort.bam.cpm']),
          ],
          tumor: [
            Float(region['X5_h3k27ac_bowtie2_hg19_sort.bam.cpm']),
            Float(region['X6_h3k27ac_bowtie2_hg19_sort.bam.cpm']),
          ],
        },
      },
      input_2percent: {
        coverage: {
          UPPP: [
            Integer(region['X3_2percent_input_bowtie2_hg19_sort.bam']),
            Integer(region['X4_2percent_input_bowtie2_hg19_sort.bam']),
          ],
          tumor: [
            Integer(region['X5_2percent_input_bowtie2_hg19_sort.bam']),
            Integer(region['X6_2percent_input_bowtie2_hg19_sort.bam']),
          ],
        },
        cpm: {
          UPPP: [
            Float(region['X3_2percent_input_bowtie2_hg19_sort.bam.cpm']),
            Float(region['X4_2percent_input_bowtie2_hg19_sort.bam.cpm']),
          ],
          tumor: [
            Float(region['X5_2percent_input_bowtie2_hg19_sort.bam.cpm']),
            Float(region['X6_2percent_input_bowtie2_hg19_sort.bam.cpm']),
          ],
        },
      },
    }
  }
end

chrom_sizes = load_chrom_sizes('source_data/hg19.chrom.sizes')

FileUtils.mkdir_p(destination_folder)
[
  {prefix_src: 'SE', prefix_dst: 'SE'},
  {prefix_src: 'E', prefix_dst: 'enhancer'},
  {prefix_src: 'P', prefix_dst: 'promoter'},
].each{|info|
  region_type = info[:prefix_dst]
  all_regions = load_regions("source_data/#{ info[:prefix_src] }_domains_EDGER_coverage_and_CPM.csv")

  [:UPPP, :tumor].each{|sample_type|
    regions = all_regions.select{|region|
      enough_coverage = region[:samples][:cpm][sample_type].any?{|cpm| cpm > 1.0 }
      # # Folder: intervals_new
      # tumor_specific = (region[:logFC] > 0)
      # uppp_specific = (region[:logFC] < 0)
      # sample_type_specific = [
      #   (sample_type == :tumor) && tumor_specific,
      #   (sample_type == :UPPP) && uppp_specific,
      # ].any?
      # sample_type_specific && region[:fdr] <= 0.1 # && enough_coverage

      # Folder: intervals_bySource
      (region[:source] == sample_type.to_s) && enough_coverage
    }

    kwargs = {sample_type: sample_type, region_type: region_type}
    store_regions_and_shades(regions, "#{destination_folder}/#{region_type}_#{sample_type}", chrom_sizes: chrom_sizes, **kwargs)
  }
}
