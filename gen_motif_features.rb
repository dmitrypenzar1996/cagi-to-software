require 'fileutils'
require 'shellwords'

def sarus_results(snvs_fn, pwm_fn, thresholds_fn, pvalue_cutoff_or_besthit:, motif_length:, snv_index:, arity:)
  package = (arity == :di) ? 'ru.autosome.di.SARUS' : 'ru.autosome.SARUS'
  cmd = [
    "cat #{snvs_fn}",
    "ruby snv2seq.rb #{snv_index} #{motif_length}",
    "java -cp sarus.jar #{package} - #{pwm_fn.shellescape} #{pvalue_cutoff_or_besthit} --pvalues-file #{thresholds_fn.shellescape} --threshold-mode pvalue --output-scoring-mode logpvalue --dont-add-flanks"
  ].join('|')
  `#{cmd}`.lines.slice_before{|l|
    l.start_with?('>')
  }.map{|chunk|
    name = chunk.first.chomp[1..-1]
    occurences = chunk.drop(1).each_with_object(Hash.new){|l, hsh|
      pvalue, pos, strand = l.chomp.split("\t")
      pvalue = Float(pvalue)
      pos = Integer(pos)
      hsh[ [pos, strand] ] = pvalue
    }
    [name, occurences]
  }.to_h
end

arity = ARGV[0].to_sym # :mono # :di
snvs_fn = ARGV[1] # 'snvs.tsv' # 'validation_snvs.tsv'
motif_features_folder = ARGV[2] # 'motif_features' # 'validation_motif_features'
ROUNDING = 5

FileUtils.mkdir_p("#{motif_features_folder}/#{arity}")
pwm_folder = "motif_collections/#{arity}/pwm/"
thresholds_folder = "motif_collections/#{arity}/thresholds/"

tfs_to_take = File.readlines('source_data/chipseq_and_motif_mat.txt').drop(1).map{|l| l.split("\t").first }

Dir.glob("#{pwm_folder}/*.{pwm,dpwm}").select{|fn|
  tfs_to_take.include?( File.basename(fn).split('.').first )
}.sort.each do |pwm_fn|
  pwm_name = File.basename(pwm_fn, File.extname(pwm_fn))
  thresholds_fn = "motif_collections/#{arity}/thresholds/#{pwm_name}.thr"
  motif_length = File.readlines(pwm_fn).size - 1
  motif_length += 1  if arity == :di

  besthits_wild_type   = sarus_results(snvs_fn, pwm_fn, thresholds_fn, motif_length: motif_length, arity: arity, snv_index: 1, pvalue_cutoff_or_besthit: 'besthit')
  all_hits_mutant_type = sarus_results(snvs_fn, pwm_fn, thresholds_fn, motif_length: motif_length, arity: arity, snv_index: 2, pvalue_cutoff_or_besthit: 1)

  File.open("#{motif_features_folder}/#{arity}/#{pwm_name}.tsv", 'w') do |fw|
    headers = ['SNV']
    headers += ["pos_1", "strand_1", "pval_1_best", "pval_2_samepos", "pos_2", "strand_2", "pval_2_best", "FC", "FC_samepos"].map{|nm| "#{pwm_name}:#{nm}" }
    fw.puts headers.join("\t")
    besthits_wild_type.each{|snv_variant_name, besthit|
      besthit_stranded_pos, besthit_logpvalue_1 = besthit.to_a.first
      pos, strand = besthit_stranded_pos
      snv_name, variant = snv_variant_name.reverse.split('%',2).map(&:reverse).reverse
      logpvalue_2_samepos = all_hits_mutant_type["#{snv_name}%2"][besthit_stranded_pos]
      (pos_2, strand_2), besthit_logpvalue_2 = all_hits_mutant_type["#{snv_name}%2"].max_by{|stranded_pos_2, logpvalue_2| logpvalue_2 }
      infos = [
        snv_name,
        pos, (strand == '+' ? 1 : -1), besthit_logpvalue_1.round(ROUNDING),
        logpvalue_2_samepos.round(ROUNDING),
        pos_2, (strand_2 == '+' ? 1 : -1), besthit_logpvalue_2.round(ROUNDING),
        (besthit_logpvalue_2 - besthit_logpvalue_1).round(ROUNDING), (logpvalue_2_samepos - besthit_logpvalue_1).round(ROUNDING)
      ]
      fw.puts infos.join("\t")
    }
  end
end
