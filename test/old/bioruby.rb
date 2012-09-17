#!/usr/local/bin/ruby
require 'bio'
 
ff = Bio::FlatFile.open(Bio::FastaFormat, ARGV[0])


qf = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])

i=0
while ((f_seq= ff.next_entry) && (q_seq = qf.next_entry))

  if f_seq.entry_id!=q_seq.entry_id
  raise "ERROR in name"
  end
  if (f_seq.seq.size!=q_seq.data.count(' ')+1)
  	raise "ERROR in sizes #{f_seq.data.size}, #{q_seq.data.count(' ')+1}"
  end

  puts f_seq.entry_id
	puts f_seq.seq
	puts q_seq.entry_id
	puts q_seq.data

	i += 1  
end


