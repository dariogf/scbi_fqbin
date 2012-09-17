#!/usr/bin/env ruby

require File.expand_path(
    File.join(File.dirname(__FILE__), %w[.. lib scbi_fastabin]))

#check args
if ARGV.count != 2
	puts "$0 fbin_file sequence_name"
	exit
end


bin_file = ARGV[0]
index_file = bin_file+'.index'
seq_name = ARGV[1]

if !File.exists?(bin_file) 
	puts "Binary file \"#{bin_file}\" doesn't exists'"
	exit
end

fb=Fastabin.new(bin_file,'r')
n,f,q=fb.read_seq(seq_name)

if n.nil?
	puts "Sequence not found"
else
	puts ">"+n
	puts f
	puts	
	puts ">"+n
	puts q
end

fb.close

