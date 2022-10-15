import DataStructures
import Shuffle
import Plots

using DataStructures
using Shuffle
using Plots; gr() #I heard somewhere using gr() is faster than standard Plots
using Random
#]add Combinatorics
#]add InspectDR
using Distributions
using Combinatorics
using InspectDR
#we want a random field Ising Model,
#First let's start with the generating all possible basis
function generator(N)
  #First instance, we want to see (1,1), (1,-1), (-1,1),(-1,-1)
  #for N=2 etc.
  a = [1,-1] #initial 2-state system for spin up or down
  states=repeat(a,N) # creates [1,-1] multiple times of N-states
  configuration = unique(collect(combinations(repeat(a,N),N)))
  #creates complete configuration for the state
#  print(configuration)
  return configuration
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
  for i=1:length(stor)
    config = getindex(stor, i) 
    #shortens the time for the for loops compared to previous
    energy = get_energy(config, h, J)
    engs[i] = energy
  end
  gs_energy = minimum(engs)
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
  #plot1 = InspectDR.Plot2D(1:2^N, energies)
  plot(1:2^N, energies)
end

#This line executes the program:
exhaustive_ising(2)


#Plot E vs configuration for all configurations
#cmd line: include("Desktop/PhysResearch/exhtIsing.jl")
