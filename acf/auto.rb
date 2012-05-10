a = Array.new
while line=gets
  a.push line.to_f
end

ave = 0.0
a.each{|v|
  ave = ave + v
}
ave = ave / a.size.to_f

var = 0.0
a.each{|v|
  var = var + (v-ave)**2
}
var = var / (a.size - 1)
std = Math.sqrt(var)
a.map!{|x| (x - ave)/std}

def acf(a,diff)
  s = a.size - diff
  r = 0.0
  for i in 0..(s-1)
    r = r+ a[i]*a[i+diff] 
  end
  r = r / s.to_f
  r
end

for i in 0..(a.size*0.01)
  puts acf(a,i)
end
