#!/usr/bin/env ruby

require 'scbi_fastabin'

# show the problem
def show_memory
 3.times { GC.start }  # try to clean up
 mem = `ps -o rss -p #{Process.pid}`[/\d+/]
 puts "Current memory:  #{mem}"
end

fb=FbinFile.new(File.expand_path('~/seqs/SRR069473.fbin'),'r',false)

fb.each do |n,s,q,e| 
  puts "@#{n}\n#{s}\n+#{n}\n#{q}"
  # gets
end

fb.close
