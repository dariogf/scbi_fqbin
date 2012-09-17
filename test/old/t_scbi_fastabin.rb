
require File.dirname(__FILE__) + '/test_helper.rb'

class TestScbiFastabin < Test::Unit::TestCase

  def setup
  	#File.delete(TEST_FILE)
  end
  
	  TEST_FILE='/tmp/fbinfile';
  	
  	SEQ_FASTA='ACTG'
		SEQ_QUAL=[25]
	  SEQ_NAME='SEQ'
	  SEQ_EXTRAS='SOME EXTRAS IN SEQ'
 
  def this_method
     caller[0] =~ /`([^']*)'/ and $1
  end

  	 
  def fill_file(n)
  	fb=Fastabin.new(TEST_FILE,'w')
  	n.times do |c|
  	  i = c+1
  		fb.add_seq(SEQ_NAME+i.to_s,SEQ_FASTA*i,(SEQ_QUAL*i*SEQ_FASTA.length).join(' '),SEQ_EXTRAS)
  	end  	

  	fb.close
  end
  
  def test_new

		if File.exists?(TEST_FILE)
 			 File.delete(TEST_FILE)
		end
		
  	fb=Fastabin.new(TEST_FILE,'w')
  	fb.close
  	
    assert(File.exists?(TEST_FILE))
   end

#  def test_add1
#		
#    # make new file and fill with data
#		fill_file(1)  	
#		
#  	fb=Fastabin.new(TEST_FILE,'r')
#    assert_equal(1,fb.count)
#    fb.close
#    
#   end
   
   def test_add100

    # make new file and fill with data
		fill_file(100)
		  	
  	fb=Fastabin.new(TEST_FILE,'r')  	  	
    assert_equal(100,fb.count)
    fb.close
 		
   end

   def test_read

    	 # make new file and fill with data
		  fill_file(100) 	
		

			fb=Fastabin.new(TEST_FILE,'r')
			
    	100.times do |c|
    		i = c+1
				n,s,q,e=fb.read_seq(SEQ_NAME+i.to_s)
				#puts n,s.length,q.split(' ').length
			  assert_equal(SEQ_NAME+i.to_s,n)
			  assert_equal(SEQ_FASTA*i,s)
			  assert_equal((SEQ_QUAL*i*SEQ_FASTA.length).join(' ') , q)
			  assert_equal(SEQ_EXTRAS , e)
			end
			 
			n,s,q,e=fb.read_seq(SEQ_NAME+'NO_EXIST')
			assert(n.nil?)
			 			 		  
		  fb.close			
   end

   
   def test_read_no_exists
   
       fill_file(100)
		   
		   fb=Fastabin.new(TEST_FILE,'r')
			 n,s,q=fb.read_seq(SEQ_NAME+'NO_EXIST')
			 assert(n.nil?)
			 
   end
   
   def test_each
   	   fill_file(100)
   	
   		 i = 1
   		 fb=Fastabin.new(TEST_FILE,'r')
			 fb.each do |n,s,q,e|
				 assert_equal(SEQ_NAME+i.to_s,n)
				 assert_equal(SEQ_FASTA*i,s)
  			 assert_equal((SEQ_QUAL*i*SEQ_FASTA.length).join(' '),q)
  			 assert_equal(SEQ_EXTRAS,e)
				 i+=1
			 end  		 
			 
			 assert_equal(i,101)
			 
  		 fb.close			
   
   end
   
   def test_each_by_index

   	   fill_file(100)

   		 i = 1
   		 fb=Fastabin.new(TEST_FILE,'r')
   		 
			 fb.each_by_index do |n,s,q,e|
				 assert_equal(SEQ_NAME+i.to_s,n)
				 assert_equal(SEQ_FASTA*i,s)
  			 assert_equal((SEQ_QUAL*i*SEQ_FASTA.length).join(' '),q)
  			 assert_equal(SEQ_EXTRAS,e)
				 i+=1
			 end  		 
			 
			 assert_equal(i,101)
  		 fb.close
   
   end
   
end
