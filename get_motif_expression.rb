require 'cgi'

motif_infos = (
  File.readlines('source_data/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv').drop(1) +
  File.readlines('source_data/HOCOMOCOv11_full_annotation_HUMAN_di.tsv').drop(1)
).map{|l|
  motif, tf, model_length, quality, rank, consensus, release, source, best_auroc_human, best_auroc_meous, peak_sets_human, peak_sets_mouse, words, family, subfamily, hgnc, entrez, uniprot_id, uniprot_ac = l.chomp.split("\t")
  {motif: motif, tf: tf, uniprot_ac: uniprot_ac, family: family, subfamily: subfamily}
}.sort_by{|infos| infos[:motif] }

samples_order = File.readlines('source_data/expression_profile.tsv').first.chomp.split("\t").drop(7).map{|s| CGI.unescape(s)[-11..-1] }
cell_lines = {
  'GBM' => ['10444-106F3', '10485-107A8'],
  'HaCaT' => ['11272-116H3'],
  'HepG2' => ['10818-111B8', '10819-111B9', '10820-111C1'],
  'HEL 92.1.7' => ['10835-111D7', '10836-111D8', '10837-111D9'],
  'HEK293T' => ['10450-106F9'],
  'K562' => ['10454-106G4', '10824-111C5', '10825-111C6', '10826-111C7'],
  'MIN6' => ['10698-109G5', '11249-116E7'],
  'SK-MEL-28' => ['10465-106H6', '10514-107E1'],
}
cell_line_sample_indices = cell_lines.map{|cell_line, samples|
  sample_indices = samples.flat_map{|sample|
    samples_order.each_with_index.select{|smp,idx|
      smp == sample # technical replicas can have the same id (e.g. 10444-106F3)
    }.map{|smp,idx| idx }.uniq
  }
  [cell_line, sample_indices]
}.to_h

transcript_expression_data = File.readlines('source_data/expression_profile.tsv').drop(3).map{|l|
  annotation, short_desc, desc, assoc_w_transcript, entrez, hgnc, uniprot_ac, *expressions = l.chomp.split("\t")
  {
    genes: short_desc.split(',').grep(/^p\d+@.+/).map{|s| s.split('@')[1] },
    uniprot_acs: uniprot_ac == 'NA' ? [] : uniprot_ac.split(',').map{|s| s.sub(/^uniprot:/,'') },
    expressions: expressions.map{|x| Float(x) },
  }
}

expression_data_by_uniprot_ac = Hash.new{|h,k| h[k] = []}
expression_data_by_gene = Hash.new{|h,k| h[k] = []}
transcript_expression_data.each_with_index{|infos, idx|
  infos[:uniprot_acs].each{|uniprot_ac|
    expression_data_by_uniprot_ac[uniprot_ac] << idx
  }

  infos[:genes].each{|gene|
    expression_data_by_gene[gene] << idx
  }
}

expression_by_gene = {}
cell_lines_order = cell_lines.keys.sort

motif_infos.each{|infos|
  tf_indices = (expression_data_by_gene[infos[:tf]] + expression_data_by_uniprot_ac[infos[:uniprot_ac]]).uniq
  tf_expression_over_samples = transcript_expression_data.values_at(*tf_indices).map{|h| h[:expressions] }.transpose.map(&:sum)
  tf_expression_over_samples = [0.0] * samples_order.size  if tf_expression_over_samples.empty?
  tf_expression_over_cell_lines = cell_lines_order.map{|cell_line|
    sample_indices = cell_line_sample_indices[cell_line]
    expressions_for_same_cell_line_samples = tf_expression_over_samples.values_at(*sample_indices)
    expressions_for_same_cell_line_samples.sum(0.0) / expressions_for_same_cell_line_samples.size
  }
  expression_by_gene[ infos[:tf] ] = {
    tf_expression_over_cell_lines: tf_expression_over_cell_lines,
    tf_expression_over_samples: tf_expression_over_samples
  }
}

File.open('motif_expressions.tsv', 'w') do |fw|
  fw.puts ['motif', *cell_lines_order, *samples_order].join("\t")
  motif_infos.each{|motif_info|
    expression_infos = expression_by_gene[ motif_info[:tf] ]
    row = [motif_info[:motif], *expression_infos[:tf_expression_over_cell_lines], *expression_infos[:tf_expression_over_samples]]
    fw.puts row.join("\t")
  }
end

motif_genes = motif_infos.map{|info| info[:tf] }.uniq
File.open('gene_expressions.tsv', 'w') do |fw|
  fw.puts ['gene', *cell_lines_order, *samples_order].join("\t")
  motif_genes.each{|gene|
    expression_infos = expression_by_gene[gene]
    row = [gene, *expression_infos[:tf_expression_over_cell_lines], *expression_infos[:tf_expression_over_samples]]
    fw.puts row.join("\t")
  }
end

gene_families = motif_infos.group_by{|info| info[:family] }.map{|family, infos| genes = infos.map{|info| info[:tf] }; [family, genes] }
File.open('motif_family_expressions.tsv', 'w') do |fw|
  fw.puts ['MotifFamily', *cell_lines_order, *samples_order].join("\t")
  gene_families.each{|family, genes|
    expressions_over_samples = genes.map{|gene| expression_by_gene[gene][:tf_expression_over_samples] }.transpose.map(&:sum)
    expressions_over_cell_lines = genes.map{|gene| expression_by_gene[gene][:tf_expression_over_cell_lines] }.transpose.map(&:sum)
    row = [family, expressions_over_cell_lines, expressions_over_samples]
    fw.puts row.join("\t")
  }
end
