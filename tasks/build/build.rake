desc "build lib"
task :build do
  chdir "libfbin_src" do
    sh "make"
  end
end

desc "Install lib"
task :install do
  chdir "libfbin_src" do
    sh "sudo make install"
  end
end


