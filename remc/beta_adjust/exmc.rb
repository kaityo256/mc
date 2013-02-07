N = 10

class Replica
  def energy(beta)
    1.0/beta
  end
end

LOOP = 100
def mcloop(b, bi,r,a)
LOOP.times{|j|
  if j%2 == 0
    (N/2).times{|i|
      i1 = bi[i*2]
      i2 = bi[i*2+1]
      b1 = b[i*2]
      b2 = b[i*2+1]
      w = (b1-b2)*(r[i1].energy(b1) - r[i2].energy(b2))
      if w > 0 || Math.exp(w) > rand
        a[i*2+1] = a[i*2+1] + 1.0/LOOP.to_f
        bi[i*2],bi[i*2+1] = bi[i*2+1],bi[i*2]
      end
    }
  else
    (N/2-1).times{|i|
      i1 = bi[i*2+1]
      i2 = bi[i*2+2]
      b1 = b[i*2+1]
      b2 = b[i*2+2]
      w = (b1-b2)*(r[i1].energy(b1) - r[i2].energy(b2))
      if w > 0 || Math.exp(w) > rand
        a[i*2+2] = a[i*2+2] + 1.0/LOOP.to_f
        bi[i*2+1],bi[i*2+2] = bi[i*2+2],bi[i*2+1]
      end
    }
  end
}
end

def beta_onestep(b,bi,r)
  a = Array.new(N){0}
  mcloop(b,bi,r,a)
  c = 0.0
  a.each{|v|
    c = c + v
  }
  c = c / (N-1).to_f
  b2 = b.clone
  (N-1).times{|i|
    b[i+1] = b[i] + (b2[i+1] - b2[i])*a[i+1]/c
  }
end

b = Array.new(N){|i| i+1}
bi = Array.new(N){|i| i}
r = Array.new(N){Replica.new()}
r.shuffle!

1000.times{|i|
  beta_onestep(b,bi,r)
}
N.times{|i|
  puts "#{i} #{b[i]}"
}
