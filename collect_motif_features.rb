class Array
  def mean(init=0.0); sum(init) / length; end
  def median; sort[length/2]; end
  def stddev; m = mean; (map{|x| (x-m) ** 2 }.sum(0.0) / length) ** 0.5; end
end

def read_expressions_table(fn)
  cell_lines = File.readlines(fn).first.chomp.split("\t").drop(1)
  groups = []
  expr_matrix = []
  File.readlines(fn).drop(1).map{|l|
    group_name, *expr_row = l.chomp.split("\t")
    groups << group_name
    expr_matrix << expr_row
  }
  [groups, cell_lines.zip(expr_matrix.transpose).to_h]
end


motif_infos = (
  File.readlines('source_data/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv').drop(1) +
  File.readlines('source_data/HOCOMOCOv11_full_annotation_HUMAN_di.tsv').drop(1)
).map{|l|
  motif, tf, model_length, quality, rank, consensus, release, source, best_auroc_human, best_auroc_meous, peak_sets_human, peak_sets_mouse, words, family, subfamily, hgnc, entrez, uniprot_id, uniprot_ac = l.chomp.split("\t")
  {motif: motif, tf: tf, uniprot_ac: uniprot_ac, family: family, subfamily: subfamily, length: model_length, rank: rank, quality: quality }
}.reject{|infos|
  infos[:quality] == 'D'
  # false
}.sort_by{|infos| infos[:motif] }

respresenative_motifs = motif_infos.group_by{|infos| infos[:family] }.map{|family, motif_infos|
  respresenative_motif = motif_infos.sort_by{|infos| [infos[:quality], infos[:rank], infos[:length]] }.first
  respresenative_motif[:motif]
}


motif_families_order, motif_family_expressions_by_cell_line = read_expressions_table('motif_family_expressions.tsv')

cell_line_by_construction = (
  File.readlines('source_data/promoters.bed').drop(1) +
  File.readlines('source_data/enhancers.bed').drop(1)
).map{|l| l.chomp.split("\t").values_at(3,6) }.to_h


cell_lines = File.readlines('motif_expressions.tsv').first.chomp.split("\t").drop(1)
expressions_by_motif = {}
File.readlines('motif_expressions.tsv').drop(1).map{|l|
  motif, *exprs = l.chomp.split("\t")
  exprs = exprs.map{|x| Float(x) }
  expressions_by_motif[motif] = cell_lines.zip(exprs).to_h
};expressions_by_motif.size

value_and_confidence_by_snv = {}
Dir.glob('source_data/train_variants/*.tsv').map{|fn|
  construction_name = File.basename(fn, '.tsv')
  File.readlines(fn).drop(1).map{|l|
    chr, pos, reference_allele, alternative_allele, value, confidence = l.chomp.split("\t")
    snv = "#{construction_name}@#{pos}@#{reference_allele}/#{alternative_allele}"
    value_and_confidence_by_snv[snv] = [value, confidence]
  }
}

ROUNDING = 3

motif_features_folder = ARGV[0] # 'motif_features' # 'validation_motif_features'

motif_features = Dir.glob("#{motif_features_folder}/*/*.tsv").sort.map{|fn|
  motif_data = File.readlines(fn).drop(1).map{|l|
    snv, pos_1, strand_1, pval_1_best, pval_2_samepos, pos_2, strand_2, pval_2_best, fc_best, fc_samepos = l.chomp.split("\t")
    datapoint = {
      pval_1_best: Float(pval_1_best).round(ROUNDING),
      pval_2_best: Float(pval_2_best).round(ROUNDING),
      pval_2_samepos: Float(pval_2_samepos).round(ROUNDING),
      fc_samepos: Float(fc_samepos).round(ROUNDING),
    }
    [snv, datapoint]
  }.to_h
  motif_name = File.basename(fn, '.tsv')
  [motif_name, motif_data]
}.reject{|motif_name, motif_data|
  quality = motif_name.split('.').last
  quality == 'D'
  # false
}.to_h

$stderr.puts 'Loaded'

motif_feature_headers = [
  'pval_1_best', 'fc_samepos',
  'pval_2_best', 'pval_2_samepos',
]

feature_names = []
features_by_snv = Hash.new{|h,k| h[k] = [] }

motif_features.each{|motif, motif_data|
  header_block = motif_feature_headers.map{|feature| "#{motif}:#{feature}" }
  feature_names.push(*header_block)
  motif_data.each{|snv, datapoint|
    construction =  snv.split('@').first
    cell_line = cell_line_by_construction[construction]

    # expr = expressions_by_motif.fetch(motif, Hash.new(0.0))[cell_line]

    fc = datapoint[:fc_samepos]
    pval = datapoint[:pval_1_best]
    pval_2_best = datapoint[:pval_2_best]
    pval_2_samepos = datapoint[:pval_2_samepos]

    infos_block = [
      pval, fc,
      pval_2_best, pval_2_samepos,
    ].map{|x| x.round(ROUNDING) }

    features_by_snv[snv].push(*infos_block)
  }
}

motifs_in_families = motif_infos.group_by{|infos|
  infos[:family]
}.map{|family, infos|
  [family, infos.map{|info| info[:motif] }]
}.sort.to_h

tbl_headers = [
  'SNV', 'construction', 'value', 'confidence',
  *feature_names,
  *motif_families_order.map{|nm| "#{nm}:family_expression"},
]

puts tbl_headers.join("\t")
features_by_snv.sort.each{|snv, features|
  # raise "Has no value&confidence for #{snv}"  unless value_and_confidence_by_snv.has_key?(snv)
  construction = snv.split('@').first
  cell_line = cell_line_by_construction[construction]
  motif_families_expressions = motif_family_expressions_by_cell_line[ cell_line ]
  
  row = [
    snv, construction, *value_and_confidence_by_snv.fetch(snv){ [nil, nil] },
    *features,
    *motif_families_expressions,
  ]
  puts row.join("\t")
}
