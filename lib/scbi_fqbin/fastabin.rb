########################################################
# Author:  Dario Guerrero & Rafael Larrosa
# 
# Encapsulates FastaBIN format routines
# 
# 
# 
########################################################

require 'zlib'
require 'scbi_fasta'

class Fastabin

	#--------------------------------------
	# CONSTANT DEFINITIONS
	#--------------------------------------

#Compression type. Valid values are Zlib::NO_COMPRESSION, Zlib::BEST_SPEED, Zlib::BEST_COMPRESSION, Zlib::DEFAULT_COMPRESSION, and an integer from 0 to 9. 
	COMPRESSION=Zlib::BEST_COMPRESSION
  SEQUENCES_PER_BLOCK=10	
  READ_BIN_REG_EXP=  /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/

  READ_REG_EXP=  /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/

	#--------------------------------------
	# Constructor
	#--------------------------------------
	def initialize(filename, mode = 'r', index_filename=nil)
		
		@filename=filename
		@index_filename = index_filename ||= filename+'.index';
		@bin_file = nil
		@need_to_regenerate_index=false

		@added_sequence_count=0
		
		# check open mode read or write
		if mode.upcase.index('W')
			bin_mode='wb'
			if File.exists?(@index_filename)
				File.delete(@index_filename)
			end
		else
		
			# if read mode, check if file exits
			if !File.exists?(filename)
				raise "File #{filename} doesn't exists'"
			end
	
			# if index doesn't exits, recreate it
			if !File.exists?(@index_filename)
			  regenerate_index
			end

			# check open mode
			bin_mode = 'rb'	
			
		end
		
		#open files
		@bin_io = File.open(filename,bin_mode)		
		@bin_file = Zlib::GzipWriter.new(@bin_io)
		
	end


	#--------------------------------------
	# Add a pair of fasta and qual files  to the Fastabin file
	#--------------------------------------
	def add_fasta_qual(fasta_file_name,qual_file_name)

			# use FastaQualFile to read fasta
				qf = FastaQualFile.new(fasta_file_name,qual_file_name)

				
				# iterate over sequences
				qf.each do	|name,fasta,qual|

          if (qf.num_seqs % 10000) == 0
						    		puts Time.now.to_s + ',' + qf.num_seqs.to_s + ':' + name
          end
#											    	
							#add them to fastabin 
							add_seq(name,fasta,qual,nil)
	
				end

				qf.close

	end
	
	#--------------------------------------
	# Add one seq to the fastabin file
	#--------------------------------------
	def add_seq(seq_name, seq_fasta, seq_qual, seq_extras=nil)

				if (@added_sequence_count % SEQUENCES_PER_BLOCK==0)
							  #@bin_file.flush
								@bin_file.close
								#@bin_io.close
								
								@bin_io = File.open(@filename,'ab')	
								@bin_file = Zlib::GzipWriter.new(@bin_io)
								puts "NEW BLOCK"
							
				end
				@added_sequence_count += 1
				
				@need_to_regenerate_index = true
				
				zipped_fasta = ''
				zipped_qual = ''
				zipped_extras = ''

				ini = @bin_file.pos
				
				# get current pos and write deflated data to fastabin format
			 	#zipped_fasta = deflate_fasta(seq_fasta.strip)
			 	zipped_fasta = seq_fasta.strip

				if !seq_qual.empty?
				  q = seq_qual.strip.split(' ')
				  q.map!{|e| (e.to_i+33).chr}
 					#zipped_qual = deflate_qual(q.join)
 					zipped_qual = q.join
 					
				  #puts q.join
				  #raise
					#zipped_qual = deflate_qual(seq_qual.strip)

				end

				if !seq_extras.nil?
					#zipped_extras = deflate_extras(seq_extras)
					zipped_extras = seq_extras
				end
		
				# write data to index file and bin (to retrieve it later if index file gets lost)
				head = "#{seq_name} #{zipped_fasta.size} #{zipped_qual.size} #{zipped_extras.size}"
				bin_index_line ="#{head.size.to_s.rjust(4)}#{head}"
				
			  #index_line ="#{seq_name} #{ini+bin_index_line.size} #{zipped_fasta.size} #{zipped_qual.size} #{zipped_extras.size}"
				#puts index_line
				#index_file.puts index_line
				
				@bin_file.write bin_index_line

#        puts zipped_fasta
#				puts zipped_qual
#				puts zipped_extras

#				puts "1F:#{zipped_fasta}@#{zipped_fasta.size.to_s}@#{zipped_fasta.length.to_s}"
#				puts "1Q:#{zipped_qual}@#{zipped_qual.size.to_s}@#{zipped_qual.length.to_s}"
#				puts "1E:#{zipped_extras}@#{zipped_extras.size.to_s}@#{zipped_extras.length.to_s}"
#				
				
				@bin_file.puts zipped_fasta
				@bin_file.puts zipped_qual
				@bin_file.puts zipped_extras
				
				
				
	end

	#--------------------------------------
	# Iterate over all sequences of a fastabin file
	#--------------------------------------	
	def each(get_fasta=true,get_qual=true,get_extras=true)
	
			#bin = ''
			@bin_file.pos=0
	
			
			while !@bin_file.eof?
			
    			head_size = @bin_file.read(4).to_i
					line = @bin_file.read(head_size)
			
					if line =~ READ_BIN_REG_EXP
			
							name = $1
							i = @bin_file.pos
							fz = $2.to_i
							qz=$3.to_i
							ez=$4.to_i
		
							@bin_file.pos = i+fz+qz+ez
							
							name,fasta,qual,extras = extract_seq(name,i,fz,qz,ez)
							
						  yield(name,fasta,qual,extras)
						
					else
						raise "Invalid index line found at each"
					end
			
			end
	
	end 
	
	#--------------------------------------
	# Iterate over all sequences of a fastabin file
	#--------------------------------------	
	def each_by_index(get_fasta=true,get_qual=true)
			  index_file = Zlib::GzipReader.open(@index_filename)
			  
				# iterate over each line of index_file				
		  	index_file.each_line do |e|

							# parse params
            if e=~ READ_REG_EXP
					
							name = $1
							i = $2.to_i
							fz = $3.to_i
							qz=$4.to_i
							ez=$5.to_i
		
							name,fasta,qual,extras = extract_seq(name,i,fz,qz,ez)
							
						  yield(name,fasta,qual,extras)
            end
							
    		end
    		index_file.close
												
	end

	#--------------------------------------
	# Read one seq from the fastabin file by name
	#--------------------------------------
	def read_seq(seq_name)

	  index_file = Zlib::GzipReader.open(@index_filename)
	  
		res = nil
		e=nil
		
		index_file.grep(/^#{seq_name}\s/) do |line|
    
				 e=line.chomp

				# parse params
				if e=~ READ_REG_EXP
					name = $1
					i = $2.to_i
					fz = $3.to_i
					qz=$4.to_i
					ez=$5.to_i

					res = extract_seq(name,i,fz,qz,ez)
				end
				break
    end
    
    index_file.close
    
	  return res 
	end
	
	
	#--------------------------------------
	# Count lines in index file. This is sequence count
	#--------------------------------------
	def count
		  index_file = Zlib::GzipReader.open(@index_filename)
    	
    	res=index_file.readlines.count

			index_file.close
			
    	return res
	end

	#--------------------------------------
	# Close files
	#--------------------------------------
	def close

		@bin_file.close if !@bin_file.closed?		
	  #@bin_io.close if !@bin_io.closed?
	  		
		regenerate_index if @need_to_regenerate_index
	end
	

private 

	#--------------------------------------
	# Extract/regenerate index from bin file
	#--------------------------------------
	
	def regenerate_index
	    #puts Time.now.to_s + " - extracting index"
	    
 			bin_io = File.open(@filename,'rb')
 			bin_file = Zlib::GzipReader.new(bin_io).to_io
		
		  index_file = Zlib::GzipWriter.open(@index_filename)
			#@bin_file.pos=0
			i=0
			while !bin_file.eof?

    			head_size = bin_file.read(4).to_i    			
					line = bin_file.read(head_size)
					
					if line =~ READ_BIN_REG_EXP
			
							name = $1
							i = bin_file.pos
							fz = $2.to_i
							qz=$3.to_i
							ez=$4.to_i
						
							index_line ="#{name} #{i} #{fz} #{qz} #{ez}"
							#puts index_line
							i += 1
							
							puts i
							index_file.puts index_line

							bin_file.pos = i+fz+qz+ez
						
					else
						raise "Invalid index line found at each"
					end
			
			end
			
			#puts Time.now.to_s + ' - done'

			index_file.close
		  bin_file.close
		  #bin_io.close
		  
	
	end
	
	def extract_seq(name,i,fz,qz,ez)
	
		#annotate current pos in files to revert changes later		
		current_bin_pos = @bin_file.pos

		res = nil
					
					# read bin data and inflate it
					@bin_file.pos = i
					
					fasta_d = @bin_file.read(fz)
#  				puts "2F:#{fasta_d}@#{fasta_d.size}"
					fasta = inflate(fasta_d)
					fasta_d=nil
			
					qual_d = @bin_file.read(qz)
					qual = inflate(qual_d)
#  				puts "2Q:#{qual_d}@#{qual_d.size}"
					qual_d=nil
					
					extras_d = @bin_file.read(ez)
					extras = inflate(extras_d)
#  				puts "2E:#{extras_d}@#{extras_d.size}"
					extras_d=nil
					
					res =[name, fasta, qual, extras]
		
				  @bin_file.pos = current_bin_pos
        
				  return res 
	end

def deflate_fasta(input)

	zipper = Zlib::Deflate.new(Zlib::BEST_COMPRESSION,15,9,Zlib::FILTERED)
	res= zipper.deflate(input, Zlib::FINISH)
	zipper.close
	
	return res
end

def deflate_qual(input)

#	zipper = Zlib::Deflate.new(Zlib::BEST_COMPRESSION,15,9,Zlib::FILTERED)
		zipper = Zlib::Deflate.new(Zlib::BEST_COMPRESSION,15,9)
	res= zipper.deflate(input, Zlib::FINISH)
	zipper.close

	return res
end

def deflate_extras(input)

	zipper = Zlib::Deflate.new(Zlib::BEST_COMPRESSION,15,9,Zlib::HUFFMAN_ONLY)
	res= zipper.deflate(input, Zlib::FINISH)
	zipper.close

	return res
end

def inflate(input)

	unzipper = Zlib::Inflate.new(15)
	res= unzipper.inflate(input)
	#unzipper.finish
	unzipper.close

	return res
end

end

