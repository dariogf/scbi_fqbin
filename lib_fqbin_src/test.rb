require 'scbi_fastq'

fqr=FastqFile.new(ARGV.shift)

r=0
ntcount=0
fqr.each do |f,q,n|
	r+=1
	ntcount+=f.length+q.length
  # puts f,q,n
end
puts "Total seqs: #{r}\n";
puts "Total NT: #{ntcount}\n";
