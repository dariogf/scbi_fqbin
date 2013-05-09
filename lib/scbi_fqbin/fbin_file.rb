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

class FbinFile
  
  CREATE_NEW_FILE=1
  APPEND_TO_FILE=2
  extend FFI::Library
  
  ffi_lib(["libfbin"])

    functions = [

      [:read_seq, [:string,:string,:pointer,:pointer,:pointer],:int],
      [:read_data_sequential, [:pointer,:pointer,:pointer,:pointer,:pointer],:int],
      [:initialize_sequential_reads,[:pointer, :string],:int],
      [:close_sequential_reads,[:pointer],:int],
      [:write_seq,[:pointer, :string, :string, :string, :string],:int],
      [:close_writes,[:pointer],:int],
      [:initialize_writes,[:pointer, :string, :int, :int, :int, :int],:int],
      [:inspect_file_data_struct,[:pointer],:void],
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
  
  
  # Initializes the file 
  def initialize(file_path, mode = 'r', qual_to_phred=true, qual_to_array=true, discretize_qual=0, flatten_qual=0, create_index=1)
    @file_path = file_path
    @open_mode=mode
    
    @qual_to_phred = qual_to_phred
    @qual_to_array = qual_to_array
    
    @to_phred = lambda{|q| q - 33}
    
    if @open_mode.index('r')
      @gzf_bin = FFI::MemoryPointer.new :pointer
      initialize_sequential_reads(@gzf_bin,@file_path)
      # inspect_file_data_struct(@gzf_bin)
      @gzf_bin = @gzf_bin.get_pointer(0)
      
    elsif @open_mode.index('w')
      @write_struct = FFI::MemoryPointer.new :pointer
      initialize_writes(@write_struct,@file_path, CREATE_NEW_FILE,discretize_qual, flatten_qual, create_index)
      @write_struct = @write_struct.get_pointer(0)
      # inspect_file_data_struct(@write_struct)
      
    elsif @open_mode.index('a')
       
       @write_struct = FFI::MemoryPointer.new :pointer
       initialize_writes(@write_struct,@file_path, APPEND_TO_FILE,discretize_qual, flatten_qual, create_index)
       @write_struct = @write_struct.get_pointer(0)
       # inspect_file_data_struct(@write_struct)
       
    else
      raise "Invalid aperture mode #{mode}. Use r/w"
    end
  end
  
  # Access a sequence by its name using the indexed read
  def read_sequence(seq_name)
    
    fastaPtr = FFIString.new
    qualPtr = FFIString.new
    extrasPtr = FFIString.new
        
    if read_seq(@file_path,seq_name,fastaPtr,qualPtr,extrasPtr)==0
      
      qual=qualPtr.to_s
      
      if @qual_to_phred
        qual=qual.each_char.map{|e| (@to_phred.call(e.ord))}
        if !@qual_to_array
          qual=qual.join(' ')
        end
      end
      
      fasta=fastaPtr.to_s
      extras = extrasPtr.to_s
      
      free_string(qualPtr)
      free_string(fastaPtr)
      free_string(extrasPtr)

      return seq_name, fasta, qual, extras
      
    else
      
      free_string(qualPtr)
      free_string(fastaPtr)
      free_string(extrasPtr)
      
      return nil
      # raise "Invalid sequence"
    end
  end
  
  # Returns the next sequence in file, or nil if no more ara available
  def next_sequence
    return get_next_seq(@gzf_bin)
  end
  
  def get_next_seq(gzf_bin)
    seq_namePtr = FFIString.new
    fastaPtr = FFIString.new
    qualPtr = FFIString.new
    extrasPtr = FFIString.new
    
    if ((r=read_data_sequential(gzf_bin, seq_namePtr, fastaPtr, qualPtr, extrasPtr))==0)
      
      qual=qualPtr.to_s
      if @qual_to_phred
        qual=qual.each_char.map{|e| (@to_phred.call(e.ord))}
        
        if !@qual_to_array
          qual=qual.join(' ')
        end
      end
      
      seq_name=seq_namePtr.to_s
      fasta=fastaPtr.to_s
      extras=extrasPtr.to_s
      
      # seq_namePtr.free
      # fastaPtr.free
      # qualPtr.free
      # extrasPtr.free
      
      free_string(seq_namePtr)
      free_string(fastaPtr)
      free_string(qualPtr)
      free_string(extrasPtr)
      
      return seq_name, fasta, qual, extras
      
    else
      free_string(seq_namePtr)
      free_string(fastaPtr)
      free_string(qualPtr)
      free_string(extrasPtr)
      return nil
    end
  end
  
  # Iterates over all sequences in file
  def each
                                                   
    gzf_bin = FFI::MemoryPointer.new :pointer
    initialize_sequential_reads(gzf_bin,@file_path)
    gzf_bin = gzf_bin.get_pointer(0)
    
     while seq=get_next_seq(gzf_bin)
       yield seq
     end
     
    close_sequential_reads(gzf_bin)
  end
  
  # Writes a sequence to file
  def write_sequence(seq_name,fasta,qual,extras)

    if qual.is_a?(Array)
        qual=qual.join(' ')
    end
      
    write_seq(@write_struct,seq_name,fasta,qual,extras)
  end
  
  
  # Closes file
  def close

    if @open_mode.index('r')
      close_sequential_reads(@gzf_bin)
    else # CREATE_NEW_FILE or APPEND_TO_FILE
      # inspect_write_struct(@write_struct)
      close_writes(@write_struct)
    end
    
  end
  
  # Returns the number of sequences in file. Iterating over it
  def count
    res=0
    each do #|n,f,q,e|
      # puts n
      res+=1
    end
    return res
  end

end
