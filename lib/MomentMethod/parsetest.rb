require 'optparse'
opt = OptionParser.new

material = 'Cu'
structure = 'jindofcc'
opt.on('-m VAL','--material','material selection, type = Cu, Au, Ag, default = Cu') {|v|
  p material = v
}
opt.on('-s VAL','--structure','structure selection, type = jindofcc, sakakifcc, default = jindofcc') {|v|
  p structure = v
}

opt.parse!(ARGV)


p material
p structure
