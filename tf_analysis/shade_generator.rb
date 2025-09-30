require 'tempfile'
require 'shellwords'

def store_regions_and_shades(regions, basename, chrom_sizes:, **kwargs)
    # should we merge shades or merging original regions is enough?
    regions_merged = merge_intervals(regions, **kwargs)
    shades_10k = regions_merged.flat_map{|region| region_shades(region, gap_size: 10_000, chrom_sizes: chrom_sizes) }
    shades_100k = regions_merged.flat_map{|region| region_shades(region, gap_size: 100_000, chrom_sizes: chrom_sizes) }
    store_regions_as_bed("#{basename}.bed", bed_sort(regions))
    store_regions_as_bed("#{basename}_shades-10k.bed", bed_sort(shades_10k))
    store_regions_as_bed("#{basename}_shades-100k.bed", bed_sort(shades_100k))
end

def merge_intervals(regions, **kwargs)
  begin
    file_sorted = Tempfile.new('sorted_regions.bed')
    file_sorted.close
    store_regions_as_bed(file_sorted.path, bed_sort(regions))
    begin
      file_merged = Tempfile.new('merged_regions.bed')
      file_merged.close
      system("bedtools merge -i #{file_sorted.path.shellescape} -c 4 -o distinct  >  #{file_merged.path.shellescape}")
      load_regions_from_bed(file_merged.path, **kwargs)
    ensure
      file_merged.unlink
    end
  ensure
    file_sorted.unlink
  end
end

def bed_sort(regions)
  regions.sort_by{|region| [region[:chr], region[:from]] }
end

def region_shades(region_info, gap_size: 10_000, chrom_sizes:)
  chr = region_info[:chr]
  from = region_info[:from]
  to = region_info[:to]
  interval_length = to - from

  left_flank = [from - gap_size - interval_length, from - gap_size]
  right_flank = [to + gap_size, to + gap_size + interval_length]
  if chrom_sizes
    left_flank = left_flank.map{|coord| [coord, 0].max }
    right_flank = right_flank.map{|coord| [coord, chrom_sizes[chr]].min }
  end

  [
    region_info.merge(from: left_flank[0],  to: left_flank[1],  id: "#{region_info[:id]}@left#{gap_size}"),
    region_info.merge(from: right_flank[0], to: right_flank[1], id: "#{region_info[:id]}@right#{gap_size}"),
  ].reject{|region|
    region[:to] <= region[:from]
  }
end

def store_regions_as_bed(filename, regions)
  File.open(filename, 'w'){|fw|
    regions.each{|info|
      row = info.values_at(:chr, :from, :to, :id) + [nil, '.']
      fw.puts(row.join("\t"))
    }
  }
end

# kwargs should be smth like {sample_type: 'UPPP', region_type: 'SE', is_specific: true}
def load_regions_from_bed(filename, **kwargs)
  File.open(filename){|f|
    f.each_line.map{|l|
      chr, from, to, id, *rest = l.chomp.split("\t")
      {chr: chr, from: Integer(from), to: Integer(to), id: id, **kwargs}
    }
  }
end

def load_chrom_sizes(filename)
  File.readlines(filename).map{|l|
    chr, sz = l.chomp.split("\t")
    [chr, sz.to_i]
  }.to_h
end
