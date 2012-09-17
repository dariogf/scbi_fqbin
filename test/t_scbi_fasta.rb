#!/usr/local/bin/ruby
require 'scbi_fasta'

ff = FastaQualFile.new(ARGV[0],ARGV[1])

ff.each do |n,f,q|
puts n
puts f
puts n
puts q
end

