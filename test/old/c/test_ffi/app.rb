require 'my_library'

o = MyLibrary.new

r=o.calculate_something(42,98.6)
puts r

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
