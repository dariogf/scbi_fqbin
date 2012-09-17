require 'ffi'

module MyModule
  extend FFI::Library
  ffi_lib "mylibrary"
  
  attach_function :calculate_something, [:int, :float], :double
  attach_function :error_code, [], :int # note empty array for functions taking zero arguments
  attach_function :create_object, [:string], :pointer
  attach_function :calculate_something_else, [:double, :pointer], :double
  attach_function :free_object, [:pointer], :void
  
end

class MyLibrary
  
  include MyModule
  
  def initialize
    puts calculate_something(1,2.3)
  end

end
