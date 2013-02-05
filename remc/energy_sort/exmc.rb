N = 10

class Replica
  def initialize(i)
    @index = i
    @energy = i*10
  end
  attr_accessor :energy, :index
end

b = Array.new(N){|i| i}
bi = Array.new(N){|i| i}
r = Array.new(N){|i| Replica.new(i)}
srand(0)
r.shuffle!

LOOP = 10
LOOP.times{|j|
  if j%2 == 0
    (N/2).times{|i|
      i1 = bi[i*2]
      i2 = bi[i*2+1]
      b1 = b[i*2]
      b2 = b[i*2+1]
      w = (b1-b2)*(r[i1].energy - r[i2].energy)
      if w > 0 || Math.exp(w) > rand
        bi[i*2],bi[i*2+1] = bi[i*2+1],bi[i*2]
      end
    }
  else
    (N/2-1).times{|i|
      i1 = bi[i*2+1]
      i2 = bi[i*2+2]
      b1 = b[i*2+1]
      b2 = b[i*2+2]
      w = (b1-b2)*(r[i1].energy - r[i2].energy)
      if w > 0 || Math.exp(w) > rand
        bi[i*2+1],bi[i*2+2] = bi[i*2+2],bi[i*2+1]
      end
    }
  end

  index = Array.new

  N.times{|i|
    index.push [i,r[bi[i]].index]
  }

  index.sort!{|i1,i2| i1[1] <=> i2[1]}
  print "#{j} "
  N.times{|i|
    #print "(#{index[i][0]},#{index[i][1]}) "
    print "#{index[i][0]} "
  }
  puts


}
