import DataStructures
import Shuffle
import Plots
using DataStructures
using Shuffle
using Plots
using Random
using Distributions
using .Threads


#Initialize system:
function initial_config(n::Int)
  config = zeros(Int64, n)
  if n%2 != 0
    println("Size n must be even.")
    return
  else
    for i=1:n
      config[i] = rand([-1, 1])
    end
  end
  #print(config)
  return config
  #we produce a randomized array of N up-or-down spins
  #e.g., N=2 may equal [1,-1]; N=4 may equal [1,-1,-1,1] as our initial configuration
end

#Create random field(s):
function gaussian_rf(N)
  nd = Normal(0, 1)
  return rand(nd, N)
end

function unit_rf(N, h)
  field = zeros(Float64, N)
  for i=1:N
    field[i] = rand([-h, h])
  end
  return field
end

#Calculate energy:
function get_energy(s, h, J)
  E0 = 0.0
  E1 = 0.0
  E2 = 0.0
  for i=1:length(s)
    #if i != length(s)
      #E0 += J*s[i]*s[i+1]
    #else
      #E0 += J*s[i]*s[1]
    #end
    for j=i+1:length(s)
       E1 += (s[i]-s[j])^2/(i-j)^2
    end
    E2 += h[i]*s[i]
  end
  E = -J*E1 - E2
  return E
end

#Add domain wall
function add_dw(config)
  n = length(config)
  if n%2 != 0
    println("size must be an even number!")
  else
    for i=Int(n/2+1):n
      config[i] = -1*config[i]
    end
  end
end

#Calculate domain wall energy cost
function dw_cost(config, h, J)
  s = copy(config)
  E_gs = get_energy(s, h, J)    # "ground state" energy
  add_dw(s)
  return get_energy(s, h, J) - E_gs
end

#Calculate magnetization  (do I need?)
function get_magnetization(config)
  M = sum(config)
  return M
end


#MC Step
function do_MC_Step(config, h, kT, J)
  N = length(config)
  E = get_energy(config, h, J)
  for i=1:N
    site = rand(1:N)
    config[site] = -1*config[site]
    #attempt to update one site of the configuration
    E_new = get_energy(config, h, J)
    #look at the hamiltonian of the configuration
    #given this updated site
    dE = E_new - E
    #look at the difference between old and new hamiltonian
    prob = exp(-dE/(kT))
    #Metropolis acceptance ratio
    r = rand(Float64)

    if min(1, prob) > r
     #if true, configuration is "updated"
      #accepted_states += 1
      E += dE
    else
      config[site] = -1*config[site]
      #if false, then configuration stays the same
    end
  end
  #M = get_magnetization(config)
  #push!(M_list, abs(M))
  #we get the magnization per spin of the updated system
end

#Metropolis Algorithm
function metropolis(config_initial, h, kT, J, mcsteps)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array

  #initial state of the configuration before a MCMC update
  #E = get_energy(config, h, J)
  # hamiltonian of a particular configuration
  Threads.@threads for i=1:mcsteps
    #For MCMC an arbitrary number steps are executed for precision
    do_MC_Step(config, h, kT, J)
    #println("MC step ", i, " at ", kT, " kT.")
  end
  return dw_cost(config, h, J)
end


#CONSTANTS:

#Number of spins in initialized configuration

#Interaction constant
const J = 1.0

#Initialize number of montecarlo steps
const mcsteps = 10000000

#temperature
const kT = 2.0

init_size = 10
iter = 20
final_size = 1000

sizes = [10, 20, 50, 100, 200, 300, 400, 500, 750, 1000]
steps = [10000, 10000, 100000, 100000, 500000, 500000, 1000000, 1000000, 2000000, 4000000]

#how should I iterate through sizes?

DW_Energy_List = Vector{Float64}()
for n in sizes
  config0 = initial_config(n)
  hg = gaussian_rf(n)
  dw_energy = metropolis(config0, hg, kT, J, popfirst!(steps))
  push!(DW_Energy_List, dw_energy)
  println("At size ", n)
end
plot(sizes, DW_Energy_List, xlabel = "System Size", ylabel = "Domain Wall Cost", linewidth = 2.5, label = "Gaussian RF")
