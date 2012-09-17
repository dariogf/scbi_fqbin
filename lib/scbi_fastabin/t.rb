#!/usr/bin/env ruby
require 'scbi_fastq'

f=FastqFile.new('/Users/dariogf/seqs/chromosomes/originals/ilumina_SRR314795.fastq')

f.each do |n,f,q,c|
  # puts n,f,q,c
end

