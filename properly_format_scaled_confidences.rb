trivariate_tag_offsets = File.readlines('source_data/trivariate_hg19.bed').map{|l|
  l.chomp.split("\t")
}.map{|chr,from,to,tag|
  [tag.upcase, Integer(from)]
}.to_h

all_lines = File.readlines('source_data/all.scaled_conf.tsv').map{|l| l.chomp.split("\t") }

header = all_lines.first
rows = all_lines.drop(1)
puts ['SNV', 'confidence_scaled', 'confidence_renormed_irf4'].join("\t")
rows.map{|row|
  tag, chr, pos, ref, alt, val, confidence, confidence_scaled, confidence_renormed_irf4 = row.values_at(0,1,2,3,4,5,6,13,15)
  if trivariate_tag_offsets.include?(tag.upcase)
    pos_relative = Integer(pos) - trivariate_tag_offsets[tag.upcase] - 1
  else
    pos_relative = pos
  end
  snv_name = "#{tag.upcase}@#{pos_relative}@#{ref}/#{alt}"
  [snv_name, confidence_scaled, confidence_renormed_irf4]
}.sort.each{|row|
  puts row.join("\t")
}