header = File.readlines('source_data/trivariate/all.with_cons.public_triv.tsv').first.chomp.split("\t")
puts ['SNV', *header.drop(10)].join("\t")
File.readlines('source_data/trivariate/all.with_cons.public_triv.tsv').drop(1).map{|l|
  #chr, pos, ref, alt, value, confidence, tag, ref_hg19, ref_matched, gerp_rs, phyloP100way, phyloP46way, phastCons46way, phast, cons100way, mz46_pre, mz100_pre, mz46_log_diff, mz100_log_diff,  mz46_log_diff_hg19, mz100_log_diff_hg19 = l.chomp.split("\t")
  chr, pos, ref, alt, value, confidence, tag, offset, ref_hg19, ref_matched, *rest = l.chomp.split("\t")
  offset = Integer(offset) -1
  snv = "#{tag.upcase}@#{offset}@#{ref}/#{alt}"
  [snv, *rest]
}.sort.each{|infos|
  puts infos.join("\t")
}