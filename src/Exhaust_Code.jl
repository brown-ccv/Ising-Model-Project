import DataStructures
import Shuffle
import Plots

using DataStructures
using Shuffle
using Plots
using Random
using Distributions

function generator(N)
  basis = [[1, -1] for i in 1:N]
  set = collect(Iterators.product(basis...))
  return vec(set)
end


function get_energy(s, h, J)
  E1 = 0.0
  E2 = 0.0
  for i=1:length(s)
    for j=i:length(s)
      if j != i
        E1 += J*(s[i]-s[j])^2/(i-j)^2
      end
    end
    E2 += h[i]*s[i]
  end
  E = E1-E2
  return E
end

function random_field(N)
  nd = Normal(0, 1)
  return rand(nd, N)
end


function get_groundstate(stor, h, J)
  groundstate = Vector{Int}()
  gs_energy = 0.0
  engs = zeros(length(stor))
  j = 1
  while !isempty(stor)
    config = pop!(stor)
    energy = get_energy(config, h, J)
    engs[j] = energy
    j+=1
    if energy <= gs_energy
      gs_energy = energy
      groundstate = config
    end
  end
  gs_data = (groundstate, gs_energy, engs)
  return gs_data
end

function exhaustive_ising(N)
  config_space = generator(N)
  #print(config_space)
  h = random_field(N)
  println("h is ", h)
  J = 0.01

  @time gs = get_groundstate(config_space, h, J)
  println()
  println("Ground state energy = ", gs[2])
  println("The ground state is ", gs[1])
  energies = gs[3]
  #print("energies are: ", energies)
  plot(1:2^N, energies)
end

#This line executes the program:
exhaustive_ising(20)


#Plot E vs configuration for all configurations
#cmd line: include("Desktop/PhysResearch/exhtIsing.jl")
