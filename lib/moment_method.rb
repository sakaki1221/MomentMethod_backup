require 'optparse'
require 'MomentMethod/version'
require 'MomentMethod/moment_method.rb'
p ARGV[0]
target_path = ARGV[0] #==nil ? './' : @argv[0]
FileUtils.cd(target_path){
  CalcMoment.new("jindofcc")
}
