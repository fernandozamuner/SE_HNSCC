Signal.trap('PIPE', 'EXIT')
require 'fileutils'
require_relative 'significance'

PRECISION = 2
REGION_TYPES = ["promoter", "enhancer", "SE"].product(["UPPP", "tumor"]).map{|region, specificity| "#{region}_#{specificity}"}
# REGION_TYPES = ["promoter", "enhancer", "SE"].product(['CellLine047', 'CellLine090', 'CellLine', 'tumor', 'UPPP']).map{|region, specificity| "#{region}_#{specificity}"}

filename = ARGV[0]
output_folder = ARGV[1]

FileUtils.mkdir_p output_folder

header, *rows = File.readlines(filename).map{|l| l.chomp.split("\t") }
bonferroni_correction = rows.size

header_part_site_notsite = [
  'odds.ratio',
  'num.sites',
  'num.not.sites',
  'neg.log10.corr.significance',
]
row_part_site_notsite = [
  :ratio,
  :a_pos,
  :a_neg,
  :corr_neglog_significance,
]

# sites vs no-sites
File.open("#{output_folder}/sites-vs-nosites.tsv", 'w'){|fw|
  output_header = ['TF']
  REGION_TYPES.each{|region_type|
    output_header << header_part_site_notsite.map{|col_name|
      "#{region_type}:site-vs-nosite:#{col_name}"
    }
    output_header << ['num.sites.shades-10k', 'num.not.sites.shades-10k', 'num.sites.shades-100k', 'num.not.sites.shades-100k'].map{|col_name|
      "#{region_type}:site-vs-nosite:#{col_name}"
    }
  }

  fw.puts output_header.join("\t")
  rows.each{|row|
    tf = row[0]
    info = [tf]
    REGION_TYPES.each{|region_type|
      significance_infos = ["shades-10k", "shades-100k"].map{|shade|
        column_names = [
          "#{region_type}",
          "#{region_type}_#{shade}",
          "neg:#{region_type}",
          "neg:#{region_type}_#{shade}",
        ]
        [shade, calc_single_tf_significance(row, header, column_names: column_names)]
      }.to_h

      significance_info = significance_infos.values.inject{|res, info_to_merge|
        merge_info(res, info_to_merge)
      }

      significance_info[:corr_neglog_significance] = correct_neglog_significance(significance_info[:neglog_significance], bonferroni_correction)
      info << significance_info.values_at(*row_part_site_notsite).map{|x| x&.round(2) }
      info << [
        *significance_infos['shades-10k'].values_at(:b_pos, :b_neg),
        *significance_infos['shades-100k'].values_at(:b_pos, :b_neg),
      ].map{|x| x&.round(2) }
    }
    fw.puts info.join("\t")
  }
}
