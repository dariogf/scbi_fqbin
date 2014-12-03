require File.dirname(__FILE__) + '/test_helper.rb'

# Ojo si va muy lento cuando se incrementa el numero de secuencias porque resulta que las longitudes de la secuencia va aumentando con cada número

class TestScbiFbinFile < Test::Unit::TestCase

  def setup
  	#File.delete(TEST_FILE)
  end
  
	  TEST_FILE='/tmp/fbinfile';
  	
    SEQ_FASTA='ACTG'
  	SEQ_FASTA400=SEQ_FASTA*100
		SEQ_QUAL='B'
		SEQ_QUAL400=(SEQ_QUAL*400*SEQ_FASTA.length)
		PHRED_QUAL=[33]
	  SEQ_NAME='SEQ'
	  SEQ_EXTRAS='SOME EXTRAS IN SEQ'
	  
	  PROGRESS=100
 
  def this_method
     caller[0] =~ /`([^']*)'/ and $1
  end

  def get_seq(i)
    # i=25
    return SEQ_FASTA*i
  end
  
  def get_qual(i)
    # i=25
    return (SEQ_QUAL*i*SEQ_FASTA.length)#.join(' ')
  end
 
  def fill_file(n)
   fb=FbinFile.new(TEST_FILE,'w',false)
   n.times do |c|
     i = c+1
     if (i%PROGRESS)==0
        puts "#{this_method}: #{i}"
     end
     
     fb.write_sequence(SEQ_NAME+i.to_s,get_seq(i),get_qual(i),SEQ_EXTRAS)
     #fb.write_sequence(SEQ_NAME+i.to_s,SEQ_FASTA400,SEQ_QUAL400,SEQ_EXTRAS)
   end
   puts "END #{this_method}"
   fb.close
  end
  
  def test_new
  
    if File.exists?(TEST_FILE)
       File.delete(TEST_FILE)
    end
      
    fb=FbinFile.new(TEST_FILE,'w',false)
    # fb.write_sequence('hola','actg','50 50 50 50','extras')
    fb.close
    
    assert(File.exists?(TEST_FILE))
  end
   
  def test_add100
    num_seqs=100
    # make new file and fill with data
    fill_file(num_seqs)
  
    fb=FbinFile.new(TEST_FILE,'r',false)
    
    assert_equal(num_seqs,fb.count)
    fb.close
  
  end
     
  def test_read_random
    num_seqs=100 #100
         # make new file and fill with data
    fill_file(100)
  
  
    fb=FbinFile.new(TEST_FILE,'r',false)
  
    num_seqs.times do |c|
      i = c+1
      n,s,q,e=fb.read_sequence(SEQ_NAME+i.to_s)
      # puts n,s,q,s.length,q.length,'_______'
      assert_equal(SEQ_NAME+i.to_s,n)
      assert_equal(get_seq(i),s)
      assert_equal(get_qual(i), q)
      assert_equal(SEQ_EXTRAS , e)
      # gets
    end
   
    fb.close      
  end
  
  def test_read_no_exists
  
     fill_file(100)
   
     fb=FbinFile.new(TEST_FILE,'r',false)
     n,s,q=fb.read_sequence(SEQ_NAME+'NO_EXIST')
     puts "H:#{n}"
     assert(n.nil?)
     fb.close
  end
  
  
   
   def test_each
       fill_file(100)
    
       i = 1
       fb=FbinFile.new(TEST_FILE,'r',false)
       
       fb.each do |n,s,q,e|
         assert_equal(SEQ_NAME+i.to_s,n)
         assert_equal(get_seq(i),s)
         assert_equal(get_qual(i),q)
         assert_equal(SEQ_EXTRAS,e)
         
         i+=1
       end       
   
       assert_equal(i,101)
   
       fb.close     
   
   end
   
  def test_each_by_index
    num_seqs=100
    
    fill_file(num_seqs)
  
    i = 1
    
    fb=FbinFile.new(TEST_FILE,'r',false)
  
    fb.each do |n,f,q,e|
      if (i%PROGRESS)==0
         puts "#{this_method}: #{i}"
      end
      assert_equal(SEQ_NAME+i.to_s,n)
      assert_equal(get_seq(i),f)
      assert_equal(get_qual(i),q)
      assert_equal(SEQ_EXTRAS,e)
      i+=1
    end
  
    assert_equal(i,num_seqs+1)
    fb.close
  
  end           
   
end
