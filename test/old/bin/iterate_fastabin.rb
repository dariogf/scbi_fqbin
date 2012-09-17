#!/usr/bin/env ruby

require File.expand_path(
    File.join(File.dirname(__FILE__), %w[.. lib scbi_fastabin]))

#check args
if ARGV.count < 2
	puts "$0 fbin_file -f|-q|-e"
	puts 
	puts "-f => Get fasta"
	puts "-q => Get qual"
	puts "-e => Get extras"	
	exit
end

bin_file = ARGV.shift
mode = ARGV.join.gsub('-','').upcase

#print mode

get_fasta=mode.index('F')
get_qual=mode.index('Q')
get_extra=mode.index('E')

index_file = bin_file+'.index'

if !File.exists?(bin_file) 
	puts "File \"#{bin_file}\" doesn't exists'"
	exit
end

# open fastabin file
fb=Fastabin.new(bin_file,'r')

# iterate over all sequences
fb.each do |n,f,q,e|
  if get_fasta 
		puts ">"+n
		puts f
	end
	
	if get_qual
		puts ">"+n
		puts q
	end
	
	if get_extra
		puts ">"+n
		puts e
	end
end

fb.close

