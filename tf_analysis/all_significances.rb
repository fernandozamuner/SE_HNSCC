Signal.trap('PIPE', 'EXIT')
require_relative 'fisher_table'

PRECISION = 2
MAX_LOGP = 1000 # replace Infinity as a log-significance with a large value (not to screw excel)
REGION_TYPES = ["promoter", "enhancer", "SE"]
# SAMPLE_TYPES = ['CellLine', 'tumor', 'UPPP', ]
SAMPLE_TYPES = ['CellLine047', 'CellLine090', 'CellLine', 'tumor', 'UPPP', ]

def merge_info(info_1, info_2)
  num_samples = [info_1[:num_samples], info_2[:num_samples]].min
  neglog_significance = [info_1[:neglog_significance], info_2[:neglog_significance]].compact.min
  a_pos = [info_1[:a_pos], info_2[:a_pos]].compact.min
  b_pos = [info_1[:b_pos], info_2[:b_pos]].compact.min
  
  if info_1[:ratio] && info_2[:ratio]
    if info_1[:ratio] <= 1 && info_2[:ratio] <= 1
      ratio = [info_1[:ratio], info_2[:ratio]].max
    elsif info_1[:ratio] >= 1 && info_2[:ratio] >= 1
      ratio = [info_1[:ratio], info_2[:ratio]].min
    else
      neglog_significance = nil
      ratio = nil
    end
  else
    ratio = nil
  end
  
  {num_samples: num_samples, ratio: ratio, neglog_significance: neglog_significance, a_pos: a_pos, b_pos: b_pos}
end

# column_names: [A+, B+, A-, B-]
def calc_single_tf_significance(row, header, column_names:)
  indices = column_names.map{|col_name| header.index(col_name) }
  a_pos, b_pos, a_neg, b_neg = row.values_at(*indices).map(&:to_i)
  fisher_table = FisherTable.by_two_classes(class_a_positive: a_pos, class_b_positive: b_pos, class_a_negative: a_neg, class_b_negative: b_neg)
  {
    neglog_significance: [-Math.log10(fisher_table.significance), 0].max,
    ratio: fisher_table.a_to_b_positive_rate_ratio,
    num_samples: fisher_table.total,
    a_pos: a_pos, b_pos: b_pos,
    a_neg: a_neg, b_neg: b_neg,
  }
end

def gen_column_names(region_type:, sample_type:, shade:)
  column_a_pos = "#{region_type}_#{sample_type}_specific"
  column_b_pos = "#{region_type}_nonspecific"

  column_a_neg = "#{region_type}_#{sample_type}_specific_#{shade}"
  column_b_neg = "#{region_type}_nonspecific_#{shade}"

  [column_a_pos, column_b_pos, column_a_neg, column_b_neg]
end

def gen_column_names_diffexpressed(region_type:, sample_type:, shade:)
  column_a_pos = "#{region_type}_tumor_specific_diffexpr-up"
  column_b_pos = "#{region_type}_tumor_specific_diffexpr-down"

  column_a_neg = "#{region_type}_tumor_specific_diffexpr-up_#{shade}"
  column_b_neg = "#{region_type}_tumor_specific_diffexpr-down_#{shade}"

  [column_a_pos, column_b_pos, column_a_neg, column_b_neg]
end

def calculate_tf_significances(rows, header, &block)
  raise 'Specify column names generator as a block'  unless block_given?
  tf_significances = {}

  REGION_TYPES.each{|region_type|
    SAMPLE_TYPES.each{|sample_type|
      rows.each{|row|
        significance_info = ["shades-10k", "shades-100k"].map{|shade|
          column_names = block.call(region_type: region_type, sample_type: sample_type, shade: shade)
          calc_single_tf_significance(row, header, column_names: column_names)
        }.inject{|res, info_to_merge|
          merge_info(res, info_to_merge)
        }
        tf = row[0]
        tf_significances[tf] ||= {}
        tf_significances[tf][region_type] ||= {}
        tf_significances[tf][region_type][sample_type] = significance_info
      }
    }
  }
  tf_significances
end

def output_significances(tf_significances, output_stream: $stdin)
  row_base = [:ratio, :num_samples, :neglog_significance]
  hdr_base = ['odds.ratio', 'total.count', 'neg.log10.corr.significance']
  hdr = REGION_TYPES.flat_map{|region_type|
    SAMPLE_TYPES.flat_map{|sample_type|
      hdr_base.map{|col_name| "#{region_type}:#{sample_type}:#{col_name}" }
    }
  }

  output_stream.puts ['TF', *hdr].join("\t")
  # tfs = tf_infos_by_region.values.flat_map(&:keys).uniq
  tf_significances.each_key{|tf|
    vals_significance = REGION_TYPES.flat_map{|region_type|
      SAMPLE_TYPES.flat_map{|sample_type|
        tf_significances[tf][region_type][sample_type].values_at(:neglog_significance)
      }
    }
    vals = REGION_TYPES.flat_map{|region_type|
      SAMPLE_TYPES.flat_map{|sample_type|
        tf_significances[tf][region_type][sample_type].values_at(*row_base)
      }
    }
    row = [tf, *vals.map{|x| x&.round(PRECISION) }]
    output_stream.puts row.join("\t") if vals_significance.any?{|x| x && x >= 2 }
  }
end

filename = ARGV[0]
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
File.open('results/all_vs_all_tumor_vs_UPPP.tsv', 'w'){|fw|
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
          "#{region_type}_UPPP",
        ]
        column_names += column_names.map{|col_name| "#{col_name}_#{shade}" }
        [shade, calc_single_tf_significance(row, header, column_names: column_names)]
      }.to_h

      significance_info = significance_infos.values.inject{|res, info_to_merge|
        merge_info(res, info_to_merge)
      }
      significance_info[:corr_neglog_significance] = [[(significance_info[:neglog_significance] || 0) - Math.log10(bonferroni_correction), 0].max, MAX_LOGP].min
      info << significance_info.values_at(*row_part).map{|x| x&.round(2) }
    }
    fw.puts info.join("\t")
  }
}


File.open('results/all_vs_all_CellLine_vs_UPPP.tsv', 'w'){|fw|
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
          "#{region_type}_UPPP",
        ]
        column_names += column_names.map{|col_name| "#{col_name}_#{shade}" }
        [shade, calc_single_tf_significance(row, header, column_names: column_names)]
      }.to_h

      significance_info = significance_infos.values.inject{|res, info_to_merge|
        merge_info(res, info_to_merge)
      }
      significance_info[:corr_neglog_significance] = [[(significance_info[:neglog_significance] || 0) - Math.log10(bonferroni_correction), 0].max, MAX_LOGP].min
      info << significance_info.values_at(*row_part).map{|x| x&.round(2) }
    }
    fw.puts info.join("\t")
  }
}


header_part_site_notsite = [
  'odds.ratio',
  'num.sites',
  'num.not.sites',
  'neg.log10.corr.significance',
]
row_part_site_notsite = [
  :ratio,
  :a_pos,
  :b_pos,
  :corr_neglog_significance,
]

# sites vs no-sites
File.open('results/sites-vs-nosites.tsv', 'w'){|fw|
  output_header = ['TF']
  REGION_TYPES.each{|region_type|
    SAMPLE_TYPES.each{|sample_type|
      output_header << header_part_site_notsite.map{|col_name|
        "#{region_type}:#{sample_type}:site-vs-nosite:#{col_name}"
      }
      output_header << ['num.sites.shades-10k', 'num.not.sites.shades-10k', 'num.sites.shades-100k', 'num.not.sites.shades-100k'].map{|col_name|
        "#{region_type}:#{sample_type}:site-vs-nosite:#{col_name}"
      }
    }
  }

  fw.puts output_header.join("\t")
  rows.each{|row|
    tf = row[0]
    info = [tf]
    REGION_TYPES.each{|region_type|
      SAMPLE_TYPES.each{|sample_type|
        significance_infos = ["shades-10k", "shades-100k"].map{|shade|
          column_names = [
            "#{region_type}_#{sample_type}",
            "neg:#{region_type}_#{sample_type}",
          ]
          column_names += column_names.map{|col_name| "#{col_name}_#{shade}" }
          [shade, calc_single_tf_significance(row, header, column_names: column_names)]
        }.to_h

        significance_info = significance_infos.values.inject{|res, info_to_merge|
          merge_info(res, info_to_merge)
        }

        significance_info[:corr_neglog_significance] = [[(significance_info[:neglog_significance] || 0) - Math.log10(bonferroni_correction), 0].max, MAX_LOGP].min
        info << significance_info.values_at(*row_part_site_notsite).map{|x| x&.round(2) }
        info << [
          *significance_infos['shades-10k'].values_at(:a_neg, :b_neg),
          *significance_infos['shades-100k'].values_at(:a_neg, :b_neg),
        ].map{|x| x&.round(2) }
      }
    }
    fw.puts info.join("\t")
  }
}
