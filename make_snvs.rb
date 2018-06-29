def parse_header(header)
  name, region = header[1..-1].split('::')
  chr, interval = region.split(':')
  from, to = interval.split('-').map(&:to_i) # [0,1)-based
  {name: name, chr: chr, from: from, to: to}
end

def parse_snv(line)
  chr, pos, ref_allele, alt_allele, *rest = line.split("\t")
  {chr: "chr#{chr}", pos: pos.to_i, reference_allele: ref_allele, alternative_allele: alt_allele}
end

def make_snv(seq, pos, alternative_allele, reference_allele: nil)
  raise "Position outside of sequence: #{pos}"  unless 0 <= pos && pos < seq.length
  raise "SNV specified reference allele `#{reference_allele}` differs from the one in reference `#{seq[pos]}` (at position #{pos}):"  unless reference_allele && seq[pos] == reference_allele
  "#{seq[0...pos]}[#{seq[pos]}/#{alternative_allele}]#{seq[(pos+1)..-1]}"
end

['source_data/promoters.fa', 'source_data/enhancers.fa'].each do |fn|
  File.readlines(fn).map(&:chomp).each_slice(2).each{|header, seq|
    infos = parse_header(header)
    seq.upcase!
    # puts header
    File.readlines("source_data/validation_variants/#{infos[:name]}.tsv").drop(1).map(&:chomp).map{|line|
      parse_snv(line)
    }.each{|snv_info|
      begin
        snv_seq = make_snv(seq, snv_info[:pos] - infos[:from] - 1, snv_info[:alternative_allele], reference_allele: snv_info[:reference_allele])
        puts "#{infos[:name]}@#{snv_info[:pos]}@#{snv_info[:reference_allele]}/#{snv_info[:alternative_allele]}\t#{'N'*30}#{snv_seq}#{'N'*30}"
      rescue => e
        $stderr.puts "Failed: #{snv_info} with #{e}"
      end
    }
  }
end
