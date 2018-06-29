trivariate_tag_offsets = File.readlines('source_data/trivariate_hg19.bed').map{|l|
  l.chomp.split("\t")
}.map{|chr,from,to,tag|
  [tag.upcase, Integer(from)]
}.to_h

# header = File.readlines('source_data/all.with_cons_upd.tsv').first.chomp.split("\t")
header = File.readlines('source_data/all.with_cons_final.tsv').first.chomp.split("\t")
puts ['SNV', *header.drop(9)].join("\t")
# File.readlines('source_data/all.with_cons_upd.tsv').drop(1).map{|l|
File.readlines('source_data/all.with_cons_final.tsv').drop(1).map{|l|
  # _chr, pos, ref, alt, _value, _confidence, tag, ref_hg19_75, ref_matched, *rest = l.chomp.split("\t")
  _chr, alt, _confidence, pos, ref, _value, tag, ref_hg19_75, ref_matched, *rest = l.chomp.split("\t")
  tag = tag.upcase
  if trivariate_tag_offsets.include?(tag)
    pos_relative = Integer(pos) - trivariate_tag_offsets[tag.upcase] - 1
  else
    pos_relative = pos
  end
  snv = "#{tag}@#{pos_relative}@#{ref}/#{alt}"
  [snv, *rest]
}.sort.each{|row|
  puts row.join("\t")
}