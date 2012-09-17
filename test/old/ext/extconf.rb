require 'mkmf'
require 'rbconfig'

dir_config("bin")

have_library("libfbin")

create_makefile("bin")
