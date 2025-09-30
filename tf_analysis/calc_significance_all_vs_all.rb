Signal.trap('PIPE', 'EXIT')
require 'fileutils'
require_relative 'significance'

PRECISION = 2
REGION_TYPES = ["promoter", "enhancer", "SE"]
# SAMPLE_TYPES = ['CellLine', 'tumor', 'UPPP', ]
SAMPLE_TYPES = ['CellLine047', 'CellLine090', 'CellLine', 'tumor', 'UPPP', ]

filename = ARGV[0]
output_folder = ARGV[1]

FileUtils.mkdir_p output_folder

header, *rows = File.readlines(filename).map{|l| l.chomp.split("\t") }
bonferroni_correction = rows.size

header_part = [
  'odds.ratio',
  'total.count',
  'neg.log10.corr.significance',
]
row_part = [
  :ratio,
  :num_samples,
  :corr_neglog_significance,
]

# all vs all
File.open("#{output_folder}/all_vs_all_tumor_vs_UPPP.tsv", 'w'){|fw|
  output_header = ['TF']
  REGION_TYPES.each{|region_type|
    output_header << header_part.map{|col_name|
      "#{region_type}:tumor-all.vs.UPPP-all:#{col_name}"
    }
  }

  fw.puts output_header.join("\t")
  rows.each{|row|
    tf = row[0]
    info = [tf]
    REGION_TYPES.each{|region_type|
      significance_infos = ["shades-10k", "shades-100k"].map{|shade|
        column_names = [
          "#{region_type}_tumor",
          "#{region_type}_tumor_#{shade}",
          "#{region_type}_UPPP",
          "#{region_type}_UPPP_#{shade}",
        ]
        [shade, calc_single_tf_significance(row, header, column_names: column_names)]
      }.to_h

      significance_info = significance_infos.values.inject{|res, info_to_merge|
        merge_info(res, info_to_merge)
      }
      significance_info[:corr_neglog_significance] = correct_neglog_significance(significance_info[:neglog_significance], bonferroni_correction)
      info << significance_info.values_at(*row_part).map{|x| x&.round(2) }
    }
    fw.puts info.join("\t")
  }
}


File.open("#{output_folder}/all_vs_all_CellLine_vs_UPPP.tsv", 'w'){|fw|
  output_header = ['TF']
  REGION_TYPES.each{|region_type|
    output_header << header_part.map{|col_name|
      "#{region_type}:CellLine-all.vs.UPPP-all:#{col_name}"
    }
  }

  fw.puts output_header.join("\t")
  rows.each{|row|
    tf = row[0]
    info = [tf]
    REGION_TYPES.each{|region_type|
      significance_infos = ["shades-10k", "shades-100k"].map{|shade|
        column_names = [
          "#{region_type}_CellLine",
          "#{region_type}_CellLine_#{shade}",
          "#{region_type}_UPPP",
          "#{region_type}_UPPP_#{shade}",
        ]
        [shade, calc_single_tf_significance(row, header, column_names: column_names)]
      }.to_h

      significance_info = significance_infos.values.inject{|res, info_to_merge|
        merge_info(res, info_to_merge)
      }
      significance_info[:corr_neglog_significance] = correct_neglog_significance(significance_info[:neglog_significance], bonferroni_correction)
      info << significance_info.values_at(*row_part).map{|x| x&.round(2) }
    }
    fw.puts info.join("\t")
  }
}
