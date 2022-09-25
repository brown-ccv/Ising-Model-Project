#import Pkg
#Pkg.add("DataStructures")
#Pkg.add("Shuffle")
#Pkg.add("Plots")
#Pkg.build("DataStructures")
#Pkg.build("Shuffle")
#Pkg.build("Plots")
using DataStructures
using Shuffle
using Plots

function generate_configs(n, stor::Stack, arr::Vector, i)

  if i == n+1
    if arr ∉ stor
      final_arr = zeros(length(arr))
      for i=1:length(arr)
        final_arr[i] = arr[i]
      end
      stor = push!(stor, final_arr)
    end
    arr = ones(n)
    return
  end

  arr[i] = 1
  generate_configs(n, stor, arr, i+1)

  arr[i] = -1
  generate_configs(n, stor, arr, i+1)
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
  if N%2 != 0
    println("N must be even!")
  else
    h = ones(N)
    for i=1:Int(N/2)
      h[i] = -1
    end
    shuffle!(h)
  end
  return h
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
  config_space = Stack{Vector{}}()
  empty_config = ones(N)
  generate_configs(N, config_space, empty_config, 1)
  #print(config_space)
  h = random_field(N)
  println("h is ", h)

  #Set J = 0.01
  gs = get_groundstate(config_space, h, 0.01)
  println()
  println("Ground state energy = ", gs[2])
  println("The ground state is ", gs[1])
  energies = gs[3]
  #print("energies are: ", energies)
  plot(1:4, energies)
end

#This line executes the program:
exhaustive_ising(2)



#Plot E vs configuration for all configurations
#cmd line: include("Desktop/PhysResearch/exhtIsing.jl")
