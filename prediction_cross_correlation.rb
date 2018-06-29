# This script used to check that predictions aren't fucked up
require 'statsample'
confidences_by_snv = File.readlines('submission.txt').drop(1).map{|l|
  chr,pos,ref,alt,construction,direction,p_direction,confidence,se,comments = l.chomp.split("\t")
  case construction
  when 'TERT-HEK293T'
    construction = 'TERT_HEK293T'
  when 'TERT-GBM'
    construction = 'TERT_GBM'
  when 'HNF4A'
    construction = 'HNF4A_P2'
  when 'MYC'
    construction = 'MYC_RS6983267'
  end
  ["#{construction}@#{pos}@#{ref}/#{alt}", Float(confidence)]
}.sort.to_h

dnases_by_snv = File.readlines('source_data/reformatted_all.with_cons.final.tsv').drop(1).map{|l|
  l.chomp.split("\t").values_at(0,5)
}.map{|snv, dnase|
  [snv, Float(dnase)]
}.to_h

ks = confidences_by_snv.keys
dnases                = dnases_by_snv.values_at(*ks).to_vector(:scale)
confidences           = confidences_by_snv.values_at(*ks).to_vector(:scale)
confidences_shuffled  = confidences_by_snv.values_at(*ks).shuffle.to_vector(:scale)
puts Statsample::Bivariate.pearson(confidences, dnases)
puts Statsample::Bivariate.pearson(confidences_shuffled, dnases)
