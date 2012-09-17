#!/usr/bin/env ruby
load './fastq_file_c.rb'

f=FastqFileC.new('/Users/dariogf/seqs/chromosomes/originals/ilumina_SRR314795.fastq')

# GC.disable

f.each do |n,f,q,c|
  # puts n,f,q,c
end

# GC.enable