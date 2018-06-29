type = Integer(ARGV[0])
motif_length = Integer(ARGV[1]) 

$stdin.each_line{|l|
  snv = l.chomp
  flank_length = motif_length - 1
  name, snv_seq = snv.split("\t")
  match = snv_seq.match(/^(?<left>[acgtn]+)\[(?<var_1>[acgtn])\/(?<var_2>[acgtn])\](?<right>[acgtn]+)$/i)
  variants = [match[:var_1], match[:var_2]]
  seq = match[:left][-flank_length .. -1] + variants[type - 1] + match[:right][0,flank_length]
  puts ">#{name}%#{type}\n#{seq}"
}