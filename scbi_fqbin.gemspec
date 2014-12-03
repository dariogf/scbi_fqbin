# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'scbi_fqbin/version'

Gem::Specification.new do |spec|
  spec.name          = "scbi_fqbin"
  spec.version       = ScbiFqbin::VERSION
  spec.authors       = ["dariogf"]
  spec.email         = ["dariogf@scbi.uma.es"]
  spec.summary       = %q{Read/write compressed fastq or fasta files in fqbin format}
  spec.description   = %q{Read/write compressed fastq or fasta files in fqbin format}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]


  spec.add_runtime_dependency "ffi"
  spec.add_development_dependency "bundler", "~> 1.7"
  spec.add_development_dependency "rake", "~> 10.0"
end
