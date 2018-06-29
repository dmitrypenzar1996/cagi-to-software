def parse_header(header)
  name, region = header[1..-1].split('::')
  chr, interval = region.split(':')
  from, to = interval.split('-').map(&:to_i) # [0,1)-based
  {name: name, chr: chr, from: from, to: to}
end

def parse_snv(line, infos)
  construction, relative_pos, alt_allele, value, confidence = line.split("\t")
  value = Float(value)
  confidence = Float(confidence)
  {chr: infos[:chr], pos: relative_pos.to_i - 1, alternative_allele: alt_allele, value: value, confidence: confidence}
end

def make_snv(seq, pos, alternative_allele, reference_allele: nil)
  raise "Position outside of sequence: #{pos}"  unless 0 <= pos && pos < seq.length
  raise "SNV specified reference allele `#{reference_allele}` differs from the one in reference `#{seq[pos]}` (at position #{pos}):"  unless !reference_allele || reference_allele && seq[pos] == reference_allele
  "#{seq[0...pos]}[#{seq[pos]}/#{alternative_allele}]#{seq[(pos+1)..-1]}"
end

['source_data/trivariate_hg19.fa'].each do |fn|
  File.readlines(fn).map(&:chomp).each_slice(2).each{|header, seq|
    infos = parse_header(header)
    seq.upcase!
    # puts header
    File.readlines("source_data/trivariate/#{infos[:name]}.tsv").map(&:chomp).map{|line|
      parse_snv(line, infos)
    }.each{|snv_info|
      begin
        snv_seq = make_snv(seq, snv_info[:pos], snv_info[:alternative_allele])
        puts "#{infos[:name]}@#{snv_info[:pos]}@#{seq[snv_info[:pos]]}/#{snv_info[:alternative_allele]}\t#{'N'*30}#{snv_seq}#{'N'*30}"
      rescue => e
        $stderr.puts "Failed: #{snv_info} with #{e}"
      end
    }
  }
end
