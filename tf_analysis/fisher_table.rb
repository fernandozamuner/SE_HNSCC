require 'rubystats'

class FisherTable
  attr_reader :class_a_positive, :class_a_negative, :class_b_positive, :class_b_negative
  def initialize(class_a_positive, class_a_negative, class_b_positive, class_b_negative)
    raise 'Class counts should be non-negative integers'  unless [class_a_positive, class_a_negative, class_b_positive, class_b_negative].all?{|x| x.is_a?(Integer) && x >= 0 }
    @class_a_positive = class_a_positive
    @class_a_negative = class_a_negative

    @class_b_positive = class_b_positive
    @class_b_negative = class_b_negative
  end

  def add_a_positive; @class_a_positive += 1; end
  def add_b_positive; @class_b_positive += 1; end

  def add_a_negative; @class_a_negative += 1; end
  def add_b_negative; @class_b_negative += 1; end

  def class_a_total
    class_a_positive + class_a_negative
  end
  def class_b_total
    class_b_positive + class_b_negative
  end

  def positive_total
    class_a_positive + class_b_positive
  end
  def negative_total
    class_a_negative + class_b_negative
  end

  def total
    class_a_total + class_b_total
  end

  def class_a_positive_rate
    (class_a_total != 0) ? (class_a_positive.to_f / class_a_total) : nil
  end
  def class_b_positive_rate
    (class_b_total != 0) ? (class_b_positive.to_f / class_b_total) : nil
  end

  def a_to_b_positive_rate_ratio
    (class_a_positive_rate && class_b_positive_rate && class_b_positive_rate != 0) ? (class_a_positive_rate / class_b_positive_rate) : nil
  end

  def a_to_b_odds_ratio
    (class_a_positive * class_b_negative).to_f / (class_a_negative * class_b_positive)  rescue nil
  end

  def b_to_a_positive_rate_ratio
    (class_a_positive_rate && class_b_positive_rate && class_a_positive_rate != 0) ? (class_b_positive_rate / class_a_positive_rate) : nil
  end

  def to_s
    "<A-: #{class_a_negative}; A+:#{class_a_positive}; B-:#{class_b_negative}; B+:#{class_b_positive}>"
  end

  def inspect; to_s; end


  def unclassified_breaks_significance?(unclassified_a: 0, unclassified_b: 0)
    # We should check whether unclassified elements can change which class has higher positive elements rate,
    # i.e. compare signs of extremal values of (rateA - rateB).
    # if (rateA_max - rateB_min) * (rateA_min - rateB_max) <= 0, than difference is insignificant.
    # So:
    # Afull = Apos + Aneg + Aundef = Atotal + Aundef
    # rateA_max = (Apos + Aundef) / Afull
    # rateA_min = (Apos) / Afull
    #
    # Not to use division we reformulate insignificance condition as:
    #   ((Apos + Aundef) * Bfull - Bpos * Afull) * (Apos * Bfull - (Bpos + Bundef) * Afull) <= 0

    a_pos = class_a_positive; a_neg = class_a_negative # just shorter aliases
    b_pos = class_b_positive; b_neg = class_b_negative
    a_total = class_a_total; a_undef = unclassified_a
    b_total = class_b_total; b_undef = unclassified_b
    a_full = a_total + a_undef
    b_full = b_total + b_undef

    ( (a_pos + a_undef) * b_full - b_pos*a_full ) * ( a_pos*b_full - (b_pos + b_undef) * a_full ) <= 0
  end

  def significance(statistical_test: Rubystats::FishersExactTest.new, unclassified_a: 0, unclassified_b: 0)
    if unclassified_a == 0 && unclassified_b == 0
      return statistical_test.calculate(class_a_positive, class_a_negative, class_b_positive, class_b_negative)[:twotail]
    end

    if unclassified_breaks_significance?(unclassified_a: unclassified_a, unclassified_b: unclassified_b)
      return 1.0 # totally insignificant difference
    end

    # # This method uses brute force to find the worst case in terms of significance.
    # # It can be very slow when there are many unclassified elements
    # (0..unclassified_a).flat_map{|a_positive_in_unclassified|
    #   a_negative_in_unclassified = unclassified_a - a_positive_in_unclassified

    #   (0..unclassified_b).map{|b_positive_in_unclassified|
    #     b_negative_in_unclassified = unclassified_b - b_positive_in_unclassified

    #     statistical_test.calculate( class_a_positive + a_positive_in_unclassified, class_a_negative + a_negative_in_unclassified,
    #                                 class_b_positive + b_positive_in_unclassified, class_b_negative + b_negative_in_unclassified
    #                               )[:twotail]
    #   }
    # }.max

    # This method is faster but can give not the exact same result because it (probably) maximizes one-tail significance, not a two-tail
    if class_a_positive * class_b_total < class_b_positive * class_a_total  # i.e.  (Apos / Atotal) < (Bpos / Btotal)
      statistical_test.calculate( class_a_positive + unclassified_a, class_a_negative,
                                  class_b_positive,                  class_b_negative + unclassified_b
                                )[:twotail]
    else
      statistical_test.calculate( class_a_positive,                  class_a_negative + unclassified_a,
                                  class_b_positive + unclassified_b, class_b_negative
                                )[:twotail]
    end
  end

  class << self
    private :new
  end
  def self.by_class_and_total(class_a_total: , class_b_total: , class_a_positive: , class_b_positive:)
    new(class_a_positive, class_a_total - class_a_positive, class_b_positive, class_b_total - class_b_positive)
  end

  def self.by_two_classes(class_a_positive: , class_b_positive: , class_a_negative: , class_b_negative:)
    new(class_a_positive, class_a_negative, class_b_positive, class_b_negative)
  end

  def self.empty
    new(0,0,0,0)
  end
end
