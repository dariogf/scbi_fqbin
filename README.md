= scbi_fqbin

* http://github.com/dariogf/scbi_fqbin

== DESCRIPTION:

scbi_fqbin is a gem to handle compressed fastq and fasta files (FQbin files) made by Rafael Larrosa, Gonzalo Claros & Dario Guerrero for the Plataforma Andaluza de Bioinformática.

FastaBin was developed to handle huge fasta and qual files in a more efficient manner. To achieve its goals, FastaBin uses some 
indexing and compression techniques in order to minimize disk accesses and thus increase the performance when reading fasta and qual files.

== FEATURES/PROBLEMS:


== SYNOPSIS:

=== Using command line tools:

* To generate a FastaBin file:

	mk_fbin fasta qual extras fbin_file

* To read a random sequence fomr FastaBin file:

	rd_seq_fbin fbin_file seqname
	
* To iterate over all sequences in FastaBin file:
	

=== Using ruby FastaBin gem:


== REQUIREMENTS:

*libfbin library

To install libfbin library, download this file http://www.scbi.uma.es/downloads/gems/libfbin.zip:

then unzip with:

unzip libfbin.zip

and install with:

make

sudo make install

By default the library will be installed to /usr/local/lib and binaries (mk_fbin, rd_seq_fbin and iterate_fbin) to /usr/local/bin

== INSTALL:

gem install scbi_fqbin


== LICENSE:

(The MIT License)

Copyright (c) 2010 Rafael Larrosa, Gonzalo Claros & Dario Guerrero

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.