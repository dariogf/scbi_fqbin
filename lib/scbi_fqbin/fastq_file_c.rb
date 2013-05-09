
# add ord method to ruby 1.8
if !String.instance_methods.include?(:ord)
   class String 
      
     def ord
       return self[0]
     end
     
   end
end
      
require 'ffi'

class FFIString<FFI::MemoryPointer
  
  def initialize
    super(:pointer,1)
    # self.write_string('a')      
    # return FFI::MemoryPointer.from_string('a'*150000000)
  end
  
  def to_s
    resPtr = self.read_pointer()

    #if null, return nil, if not return string
    return resPtr.null? ? nil : resPtr.read_string()
  end
  
  def inspect
     self.to_s
  end
  
end



class FastqFileC

  attr_accessor :num_seqs


  extend FFI::Library
  
  ffi_lib(["libfbin"])

    functions = [

      [:get_next_seq_fastq, [:pointer,:pointer,:pointer,:pointer,:pointer],:int],
      [:open_file,[:string, :pointer],:int],
      [:close_file,[:pointer],:int],
      [:free_string,[:pointer],:int]
      
      
    ]

    functions.each do |func|
      begin
        attach_function(*func)
        private func[0]
      rescue Object => e
        puts "Could not attach #{func}, #{e.message}"
      end
    end


    def open_fastq()
       @fastq_file = FFI::MemoryPointer.new :pointer
       # puts @fastq_file.address
        open_file(@fasta_file_name,@fastq_file)
        # inspect_file_data_struct(@gzf_bin)
        # puts @fastq_file.address
        
        @fastq_file = @fastq_file.get_pointer(0)
        # puts "2"
        # puts @fastq_file.address
        # if @fastq_file.null?
        #   puts "ES NULLLLLLL"
        # end
        
    end

  #------------------------------------
  # Initialize instance
  #------------------------------------
  def initialize(fasta_file_name, mode='r', fastq_type = :sanger, qual_to_array=true, qual_to_phred=true)

  	@fasta_file_name=fasta_file_name
  	
    if mode.upcase.index('W')
      @fastq_file = File.open(fasta_file_name,'w')
    elsif mode.upcase.index('A')
      if !File.exist?(fasta_file_name)
      	raise "File #{fasta_file_name} doesn't exists" 
      end
    	
      @fastq_file = File.open(fasta_file_name,'a')
    else #read only
      if !File.exist?(fasta_file_name)
      	raise "File #{fasta_file_name} doesn't exists" 
      end
    	
      if fasta_file_name.is_a?(IO)
    	  @fastq_file = fasta_file_name
      else
        # @fastq_file = File.open(fasta_file_name,'r')
        @namePtr = FFIString.new
        @fastaPtr = FFIString.new
        @qualPtr = FFIString.new
        @extrasPtr = FFIString.new
        
        open_fastq
      end
    end
    
    @mode = mode
    @num_seqs = 0
    @fastq_type=fastq_type
    
    #  S - Sanger        Phred+33,  raw reads typically (0, 40)
    #  X - Solexa        Solexa+64, raw reads typically (-5, 40)
    #  I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    #  J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    # > >>> def solexa_quality_from_phred(phred_quality) :
    # > ...     return 10*log(10**(phred_quality/10.0) - 1, 10)
    # > ...
    # > >>> solexa_quality_from_phred(90)
    # > 89.999999995657035
    # > >>> solexa_quality_from_phred(50)
    # > 49.99995657033466
    # > >>> solexa_quality_from_phred(10)
    # > 9.5424250943932485
    # > >>> solexa_quality_from_phred(1)
    # > -5.8682532438011537
    # > >>> solexa_quality_from_phred(0.1)
    # > -16.32774717238372
    # > 
    # > >>> def phred_quality_from_solexa(solexa_quality) :
    # > ...     return 10*log(10**(solexa_quality/10.0) + 1, 10)
    # > ...
    # > >>> phred_quality_from_solexa(90)
    # > 90.000000004342922
    # > >>> phred_quality_from_solexa(10)
    # > 10.41392685158225
    # > >>> phred_quality_from_solexa(0)
    # > 3.0102999566398116
    # > >>> phred_quality_from_solexa(-20)
    # > 0.043213737826425784
    
    
    #sanger by default
    @to_phred = lambda{|q| q - 33}
    @from_phred = lambda{|q| (q+33).chr}
    
    if @fastq_type == :ilumina
        @to_phred = lambda{|q| q - 64}
        # @from_phred = lambda{|q| (q+64).chr}
        
    elsif @fastq_type == :solexa
       # 
       # solexa to phred quals
       
       @to_phred = lambda{|q| (10*Math.log(10**(q/10.0)+1,10)).round}
       # @from_phred = lambda{|q| (10*Math.log(10**(q/10.0)-1,10)).round.chr}
       
       #phred to solexa quals
       
    end
    
    @qual_to_array = qual_to_array
    
    @qual_to_phred = qual_to_phred
    
  end
  
  def close
    free_string(@namePtr)
    free_string(@qualPtr)
    free_string(@fastaPtr)
    free_string(@extrasPtr)
    
    close_file(@fastq_file)
  end
 
  
  #------------------------------------
  # Iterate over all sequences
  #------------------------------------
  def each
        
    rewind

	  n,f,q,c=next_seq

    while (!n.nil?)
			yield(n,f,q,c)
			n,f,q,c=next_seq
    end

  	rewind
  	
  end

  # goto first position in file
  def rewind
     
     @num_seqs = 0;
     close_file(@fastq_file)
     open_fastq
     # @fastq_file.pos=0
    
  end

  #------------------------------------
  # Get next sequence
  #------------------------------------
  def next_seq
    #init variables
    
    # namePtr = FFIString.new
    # fastaPtr = FFIString.new
    # qualPtr = FFIString.new
    # extrasPtr = FFIString.new
    
      
    if get_next_seq_fastq(@fastq_file,@namePtr,@fastaPtr,@qualPtr,@extrasPtr)==1
      
      seq_name=@namePtr.to_s
      qual=@qualPtr.to_s
      
      # if @qual_to_phred
      #   qual=qual.each_char.map{|e| (@to_phred.call(e.ord))}
      #   if !@qual_to_array
      #     qual=qual.join(' ')
      #   end
      # end
      
      fasta=@fastaPtr.to_s
      extras = @extrasPtr.to_s
      
      # free_string(namePtr)
      # free_string(qualPtr)
      # free_string(fastaPtr)
      # free_string(extrasPtr)
      
      return seq_name, fasta, qual, extras
      
    else

      # free_string(namePtr)
      # free_string(qualPtr)
      # free_string(fastaPtr)
      # free_string(extrasPtr)
      
      return nil
      # raise "Invalid sequence"
    end
    # res = read_fastq
    # return res
  end
  
  # write sequence to file in sanger format
  def write_seq(seq_name,seq_fasta,seq_qual,comments='')
	  name = ""
	  
		@fastq_file.puts("@#{seq_name} #{comments}")
		@fastq_file.puts(seq_fasta)
		@fastq_file.puts("+#{seq_name} #{comments}")
		
		if seq_qual.is_a?(Array)
		  @fastq_file.puts(seq_qual.map{|e| @from_phred.call(e)}.join)
	  else
		  @fastq_file.puts(seq_qual.split(/\s+/).map{|e| @from_phred.call(e.to_i)}.join)
		end
    
  end

  
  # creates fastq otuput in sanger format
  def self.to_fastq(seq_name,seq_fasta,seq_qual,comments='')
    
    res=[]
    
	  name = ""
	  
		res << ("@#{seq_name} #{comments}")
		res << (seq_fasta)
		res << ("+#{seq_name} #{comments}")
		
		if @qual_to_phred
  		if seq_qual.is_a?(Array)
  		  res<<(seq_qual.map{|e| (e+33).chr}.join)
  	  else
  		  res<<(seq_qual.split(/\s+/).map{|e| (e.to_i+33).chr}.join)
  		end
		else
		  res << seq_qual
	  end
	  
    return res
  end
  
  def with_qual?
    true
  end
  
  
  private 
  
  #------------------------------------
  #  Read one sequence in fastq
  #------------------------------------
  # @GEM-108-D02
  # AAAAGCTGG
  # +
  # :::::::::

  def read_fastq

    seq_name = nil
    seq_fasta = nil
    seq_qual = nil
    comments = nil
    
    reading = :fasta
    
    if !@fastq_file.eof
      
      begin
        #read four lines
        name_line = @fastq_file.readline.chomp
        seq_fasta = @fastq_file.readline.chomp
        name2_line = @fastq_file.readline.chomp
        seq_qual = @fastq_file.readline.chomp

      
        # parse name
        if name_line =~ /^@\s*([^\s]+)\s*(.*)$/
          # remove comments
          seq_name = $1
          comments=$2
        else
          raise "Invalid sequence name in #{name_line}"
        end
      
        # parse fasta
        seq_fasta.strip! if !seq_fasta.empty?

        # parse qual_name
    
        if !seq_name.nil? && !seq_qual.empty?

           @num_seqs += 1
       
           if @qual_to_phred
             seq_qual=seq_qual.each_char.map{|e| (@to_phred.call(e.ord))}

             if !@qual_to_array
                 seq_qual=seq_qual.join(' ')
             end
           end
       
        end
      rescue EOFError
        raise "Bad format in FastQ file"
      end
    end
    
    return [seq_name,seq_fasta,seq_qual,comments]
  end
  
  
end
