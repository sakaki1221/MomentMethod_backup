require "bundler/gem_tasks"
require "rspec/core/rake_task"
require 'yard'

task :default do
  system 'rake -T'
end

desc "make documents by yard"
task :yard do
  files = Dir.entries('docs')
  files.each{|file|
    name=file.split('.')
    if name[1]=='hiki' then
      p command="hiki2md docs/#{name[0]}.hiki > MomentMethod.wiki/#{name[0]}.md"
      system command
    end
  }
  system "cp MomentMethod.wiki/readme_ja.md README.md"
  system "cp MomentMethod.wiki/readme_ja.md MomentMethod.wiki/Home.md"
  system "cp docs/*.gif MomentMethod.wiki"
  system "cp docs/*.gif doc"
  YARD::Rake::YardocTask.new
end

desc "run spec for all members."
task :spec do
         RSpec::Core::RakeTask.new(:spec)
end

desc "setenv for release from Kwansei gakuin."
task :setenv do
  p command='setenv HTTP_PROXY http://proxy.ksc.kwansei.ac.jp:8080'
  system command
  p command='setenv HTTPS_PROXY http://proxy.ksc.kwansei.ac.jp:8080'
  system command
end
