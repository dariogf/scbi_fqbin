require File.dirname(__FILE__) + '/test_helper.rb'

class TestFillfile < Test::Unit::TestCase

  def setup
  end
  
	TEST_FILE='/tmp/fbinfile';
  	
  	SEQ_FASTA='ACTG'
		SEQ_QUAL=[25]
	  SEQ_NAME='SEQ'
	  SEQ_EXTRAS='SOME EXTRAS IN SEQ'
  	 
  	 
  def fill_file(n)
  	fb=Fastabin.new(TEST_FILE,'w')
  	n.times do |c|
  	  i = c+1
  		fb.add_seq(SEQ_NAME+i.to_s,SEQ_FASTA*i,(SEQ_QUAL*i*SEQ_FASTA.length).join(' '),SEQ_EXTRAS)  		
  	end

  	fb.close
  end
  
#
#  def test_add100
#		
#    # make new file and fill with data
#		fill_file(100)  	
#		
##  	fb=Fastabin.new(TEST_FILE,'r')
##    assert(fb.count == 100)
##    fb.close
#  	
#   end
   
     def test_read

    	 # make new file and fill with data
		  fill_file(10) 	
		

			fb=Fastabin.new(TEST_FILE,'r')
			
    	10.times do |c|
    		i = c+1
				n,s,q,e=fb.read_seq(SEQ_NAME+i.to_s)
				#puts n,s.length,q.split(' ').length
			  assert(SEQ_NAME+i.to_s==n)
			  assert(SEQ_FASTA*i==s)
			  assert((SEQ_QUAL*i*SEQ_FASTA.length).join(' ')==q)
			  assert(SEQ_EXTRAS==e)
			end
			 
			n,s,q,e=fb.read_seq(SEQ_NAME+'NO_EXIST')
			assert(n.nil?)
			 			 		  
		  fb.close			
   end

   
   
   
   
end
