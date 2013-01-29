def v(x)
  return x**2-x**3
end

def mc(beta)
  xstep = 0.2
  x = 0.0
  step = 0
  while(true)
    nx = x + (rand*2-1) * xstep
    next if nx < 0.0 
    de = v(nx) - v(x)
    if de < 0 || Math.exp(-de*beta) > rand
     x = nx
    end
    if (x > 1)
      return step
    end
    step = step + 1
  end
end

AVE = 10
AVE2 = 100
def mc_beta(beta)
sum = 0.0
var = 0.0
AVE.times{
  a = 0.0
  AVE2.times{
    a = a +  mc(beta)
  }
  a = a / AVE2.to_f
  sum = sum + a
  var = var + a * a
}
sum = sum / AVE.to_f
var = var / AVE.to_f
var = var - sum**2
puts "#{beta} #{sum} #{var**0.5}"
end

bs = 1.0
be = 50.0
ND = 10
(ND+1).times{|i|
  beta = (be-bs)/ND.to_f*i + bs
  mc_beta(beta)
}

