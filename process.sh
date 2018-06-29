cat source_data/promoters.bed | tail -n+2 | cut -f1-4 | bedtools slop -b 0 -g source_data/hg19.genome | bedtools getfasta -fi source_data/hg19.fa -bed - -fo source_data/promoters.fa -name+
cat source_data/enhancers.bed | tail -n+2 | cut -f1-4 | bedtools slop -b 0 -g source_data/hg19.genome | bedtools getfasta -fi source_data/hg19.fa -bed - -fo source_data/enhancers.fa -name+
cat source_data/trivariate_hg19.bed | bedtools slop -b 0 -g source_data/hg19.genome | bedtools getfasta -fi source_data/hg19.fa -bed - -fo source_data/trivariate_hg19.fa -name+
cat source_data/trivariate_mm9.bed | bedtools slop -b 0 -g source_data/mm9.genome  | bedtools getfasta -fi source_data/mm9.fa  -bed - -fo source_data/trivariate_mm9.fa -name+
ruby make_snvs.rb > snvs.tsv
ruby make_trivariate_snvs.rb > trivariate_snvs.tsv
java -cp ape.jar ru.autosome.perfectosape.SNPScan motif_collections/mono/pwm/ snvs.tsv  --precalc motif_collections/mono/thresholds/ --log-fold-change > perfectos_results.tsv
java -cp ape.jar ru.autosome.perfectosape.di.SNPScan motif_collections/di/pwm/ snvs.tsv  --precalc motif_collections/di/thresholds/ --log-fold-change > perfectos_results_di.tsv

wget -O source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz  http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz
( \
  seq 1 7; \
  zcat source_data/robust_phase1_pls_2.tpm.desc121113.osc.txt.gz | grep -Pe '^##ColumnVariables' | cat -n \
    | grep -Pe '10818-111B8|10819-111B9|10820-111C1|10835-111D7|10836-111D8|10837-111D9|10450-106F9|10454-106G4|10824-111C5|10825-111C6|10826-111C7|10444-106F3|10485-107A8|10465-106H6|10514-107E1|11272-116H3|11249-116E7|10698-109G5' \
    | sed -re 's/^\s*([0-9]+)\s+/\1\t/' | cut -f1 \
) \
  | tr '\n' ',' \
  | xargs -I{} echo 'zcat source_data/robust_phase1_pls_2.tpm.desc121113.osc.txt.gz | grep -vPe "^#" | cut -f {}' \
  | sed -re 's/,$//' \
  | bash \
  > source_data/expression_profile.tsv

# ruby reformat_conservative_scores.rb > source_data/reformatted_all.with_cons_cage_dnase.tsv
# ruby reformat_conservative_scores_trivariate.rb > source_data/trivariate/reformatted_all.with_cons.public_triv.tsv
ruby reformat_conservative_scores.rb > source_data/reformatted_all.with_cons.final.tsv

ruby gen_motif_features.rb  mono  snvs.tsv  motif_features
ruby gen_motif_features.rb  di  snvs.tsv  motif_features
ruby gen_motif_features.rb  mono  validation_snvs.tsv  validation_motif_features
ruby gen_motif_features.rb  di  validation_snvs.tsv  validation_motif_features

ruby get_motif_expression.rb

ruby collect_motif_features.rb  motif_features > features/motif_unnormed_expression_families.tsv
ruby collect_motif_features.rb  validation_motif_features > features/validation_motif_unnormed_expression_families.tsv
ruby collect_motif_features_trivariate.rb > features/motif_unnormed_features_trivariate.tsv


sed -ire 's/MYC_rs6983267/MYC_RS6983267/' features/validation_motif_unnormed_expression_families.tsv
sed -ire 's/MYC_rs6983267/MYC_RS6983267/' features/motif_unnormed_expression_families.tsv

ruby properly_format_scaled_confidences.rb > scaled_confidences.tsv
cat source_data/all.scaled_conf.tsv | cut -f 1,8-12 | uniq | ruby -e 'puts readline; readlines.each{|x| puts x.upcase }' > source_data/scaling.tsv

# join --header -t $'\t' -j 1 \
#    <( head -1 scaled_confidences.tsv; cat scaled_confidences.tsv | tail -n+2 | sort -k1,1 ) \
#    <( head -1 features/motif_unnormed_expression_families.tsv; cat features/motif_unnormed_expression_families.tsv | tail -n+2 | sort -k1,1 ) \
#    > features/motif_unnormed_expression_families_scaled_confidence.tsv
# join --header <( join --header features/motif_unnormed_expression_families_scaled_confidence.tsv source_data/reformatted_all.with_cons.upd.tsv -t $'\t' ) train_val.blocks -t $'\t' > features/motif_unnormed_expression_families_conservativity_cage_dnase_blocks_scaled_confidence.tsv
join --header <( join --header features/motif_unnormed_expression_families.tsv source_data/reformatted_all.with_cons.final.tsv -t $'\t' ) train_val.blocks -t $'\t' > features/motif_unnormed_expression_families_conservativity_cage_dnase_blocks.tsv
# join --header <( join --header features/validation_motif_unnormed_expression_families.tsv source_data/reformatted_all.with_cons.final.tsv -t $'\t' ) train_val.blocks -t $'\t' > features/validation_motif_unnormed_expression_families_conservativity_cage_dnase_blocks.tsv
join --header features/validation_motif_unnormed_expression_families.tsv source_data/reformatted_all.with_cons.final.tsv -t $'\t' > features/validation_motif_unnormed_expression_families_conservativity_cage_dnase.tsv

# join --header -t $'\t' -j 1 \
#    <( head -1 scaled_confidences.tsv; cat scaled_confidences.tsv | tail -n+2 | sort -k1,1 ) \
#    <( head -1 features/motif_unnormed_features_trivariate.tsv; cat features/motif_unnormed_features_trivariate.tsv | tail -n+2 | sort -k1,1 ) \
#    > features/motif_unnormed_features_scaled_confidence_trivariate.tsv
# join --header features/motif_unnormed_features_scaled_confidence_trivariate.tsv source_data/reformatted_all.with_cons.upd.tsv -t $'\t' > features/motif_unnormed_features_trivariate_with_cons_cage_dnase_scaled_confidence.tsv
join --header features/motif_unnormed_features_trivariate.tsv source_data/reformatted_all.with_cons.final.tsv -t $'\t' > features/motif_unnormed_features_trivariate_with_cons_cage_dnase.tsv



find motif_collections/mono/pwm/ -xtype f -iname '*.pwm' | xargs -n1 basename -s .pwm | xargs -n1 -I{MOTIF} echo "java -cp sarus.jar ru.autosome.SARUS source_data/promoters.fa motif_collections/mono/pwm/{MOTIF}.pwm 1.0 --pvalues-file motif_collections/mono/thresholds/{MOTIF}.thr --threshold-mode pvalue --output-scoring-mode pvalue | ruby -e 'cutoffs = [0.1,0.01,0.005,0.002,0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, 0.000005, 0.000002, 0.000001]; puts [\"construction\", *cutoffs.map{|pv| \"cutoff:#{pv}\" }].join(\"\t\"); \$stdin.readlines.slice_before(/^>/).map{|hdr, *rest| pvals = rest.map{|site| Float(site.chomp.split(\"\t\").first) }; construction = hdr[1..-1].split(\"::\").first; counts = cutoffs.map{|cutoff| pvals.count{|pv| pv <= cutoff } }; puts [construction, *counts].join(\"\t\") }' > motif_occurences/{MOTIF}.tsv"

find motif_collections/di/pwm/ -xtype f -iname '*.dpwm' | xargs -n1 basename -s .dpwm | xargs -n1 -I{MOTIF} echo "java -cp sarus.jar ru.autosome.di.SARUS source_data/promoters.fa motif_collections/di/pwm/{MOTIF}.dpwm 1.0 --pvalues-file motif_collections/di/thresholds/{MOTIF}.thr --threshold-mode pvalue --output-scoring-mode pvalue | ruby -e 'cutoffs = [0.1,0.01,0.005,0.002,0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, 0.000005, 0.000002, 0.000001]; puts [\"construction\", *cutoffs.map{|pv| \"cutoff:#{pv}\" }].join(\"\t\"); \$stdin.readlines.slice_before(/^>/).map{|hdr, *rest| pvals = rest.map{|site| Float(site.chomp.split(\"\t\").first) }; construction = hdr[1..-1].split(\"::\").first; counts = cutoffs.map{|cutoff| pvals.count{|pv| pv <= cutoff } }; puts [construction, *counts].join(\"\t\") }' > motif_occurences/{MOTIF}.tsv"
