#!/usr/bin/env ruby

require 'scbi_fasta'

			# use FastaQualFile to read fasta
				qf = FastaQualFile.new(ARGV[0],ARGV[1])

				
				# iterate over sequences
				qf.each do	|name,fasta,qual| 				

					    	puts "> #{name}"
					    	puts fasta
					    	puts "> #{name}"
					    	puts qual
					    	
				end

				qf.close

