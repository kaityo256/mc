N = 10
LOOP = 100
$P = 0.6 # Desired Exchange Ratio

def energy(b)
  return 1.0/b
end
def mcloop(b,a)
LOOP.times{|j|
    (N/2).times{|i|
      b1 = b[i*2]
      b2 = b[i*2+1]
      w = (b1-b2)*(energy(b1) - energy(b2))
      if w > 0 || Math.exp(w) > rand
        a[i*2+1] = a[i*2+1] + 1.0/LOOP.to_f
      end
    }
    (N/2-1).times{|i|
      b1 = b[i*2+1]
      b2 = b[i*2+2]
      w = (b1-b2)*(energy(b1) - energy(b2))
      if w > 0 || Math.exp(w) > rand
        a[i*2+2] = a[i*2+2] + 1.0/LOOP.to_f
      end
    }
}
end

def beta_onestep(b)
  a = Array.new(N){0}
  mcloop(b,a)
  b2 = b.clone
  (N-1).times{|i|
    b[i+1] = b[i] + (b2[i+1] - b2[i])*a[i+1]/$P
  }
end

b = Array.new(N){|i| i+1}

1000.times{|i|
  beta_onestep(b)
}
N.times{|i|
  puts "#{i} #{b[i]}"
}
