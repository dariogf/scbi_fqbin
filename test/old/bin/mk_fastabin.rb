#!/usr/bin/env ruby

require File.expand_path(File.join(File.dirname(__FILE__), %w[.. lib scbi_fqbin]))

require 'zlib'

if ARGV.count != 3
	puts "$0 fasta_file qual_file out_file"
	exit
end

fasta_file = ARGV[0]
qual_file = ARGV[1]
output_name = ARGV[2] ||= File.basename(fasta_file,File.extname(fasta_file))+'.fbin'

fb=Fastabin.new(output_name,'wb')

fb.add_fasta_qual(fasta_file,qual_file)

fb.close


