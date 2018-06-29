ruby make_snvs.rb > snvs.tsv

ruby gen_motif_features.rb  mono  snvs.tsv  motif_features
ruby gen_motif_features.rb  di  snvs.tsv  motif_features

ruby collect_motif_features.rb  motif_features > features/motif_unnormed_expression_families.tsv