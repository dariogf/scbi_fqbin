require 'lib/scbi_fastabin/fbin_file'

file='c/sample/f1.fbin'

o = FbinFile.new(file,'r')


name,fasta,qual,extras=o.next_sequence
while !name.nil?
  puts "Name:#{name}, fasta: #{fasta}"
  name,fasta,qual,extras=o.next_sequence
end

puts "="*20

o.each do |name,fasta,qual,extras|
  puts "Name:#{name}, fasta: #{fasta}"
end

puts "="*20

seq='FX9YN3P05C43XJ'
name,fasta,qual,extras=o.read_sequence(seq)
puts "Name:#{name}, fasta: #{fasta}"

o.close


# r=o.calculate_something(42,98.6)
# puts r

# c = MyLibrary.calculate_something(42, 98.6) # note FFI handles literals just fine
# 
# if ( (errcode = MyLibrary.error_code()) != 0)
#   puts "error calculating something: #{errcode}"
#   exit 1
# end
# 
# objptr = MyLibrary.create_object("my object") # note FFI handles string literals as well
# d = MyLibrary.calculate_something_else(c, objptr)
# MyLibrary.free_object(objptr)
# 
# puts "calculated #{d}"
