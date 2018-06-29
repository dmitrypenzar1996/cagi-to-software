$all_mono_motifs = Dir.glob('motif_collections/mono/pwm/*.pwm').map{|fn| File.basename(fn, '.pwm') }
$mono_motifs_by_uniprot = $all_mono_motifs.group_by{|motif| motif.split('.').first }
$all_di_motifs = Dir.glob('motif_collections/di/pwm/*.dpwm').map{|fn| File.basename(fn, '.dpwm') }
$mono_corresponding_to_di = $all_di_motifs.map{|di_motif|
  uniprot = di_motif.split('.').first
  mono_motif = $mono_motifs_by_uniprot[uniprot].detect{|motif| motif.split('.')[2] == '0' }
  [di_motif, mono_motif]
}.to_h

def take_equivalent_mono_motif(motif_name)
  if motif_name.match(/\.H11MO\./)
    motif_name
  else
    $mono_corresponding_to_di[motif_name]
  end
end

def symmetrize_distances!(matrix)
  matrix.each_key{|k|
    matrix[k][k] = 0.0
  }
  matrix.each_key{|k1|
    matrix[k1].each_key{|k2|
      matrix[k1][k2] = matrix[k2][k1]  unless matrix[k1][k2]
    }
  }
  matrix
end

def read_matrix(fn)
  lns = File.readlines(fn).map{|l| l.chomp.split("\t") }
  motifs_1 = lns.first.drop(1)
  matrix = lns.drop(1).map{|row|
    motif_2, *rest = row
    vals = rest.map{|x| x == 'x' ? nil : Float(x) }
    [motif_2, motifs_1.zip(vals).to_h]
  }.to_h
  symmetrize_distances!(matrix)
end

def distance_to_similarity!(matrix)
  matrix.each_key{|k1|
    matrix[k1].each_key{|k2|
      matrix[k1][k2] = 1.0 - matrix[k1][k2]
    }
  }
  matrix
end

$similarity_matrix = distance_to_similarity!(read_matrix('HOCOMOCOv11_full_distance_matrix_HUMAN_mono.txt'))

def motif_is_around?(motif_list, motif, similarity_threshold)
  motif_mono_equivalent = take_equivalent_mono_motif(motif)
  motif_list.map{|already_choosen_motif|
    $similarity_matrix[motif_mono_equivalent][take_equivalent_mono_motif(already_choosen_motif)]
  }.any?{|x|
    x >= similarity_threshold
  }
end

similarity_threshold = Float(ARGV[0] || 0.05)

motif_list = $stdin.readlines.map(&:strip).reject(&:empty?)
result = []
motif_list.each{|motif|
  unless motif_is_around?(result, motif, similarity_threshold)
    result << motif
    puts motif
  end
}
