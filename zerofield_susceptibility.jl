import DataStructures
import Shuffle
import Plots
using DataStructures
using Shuffle
using Plots
using Random
using Distributions
using .Threads


#takes in an integer n
#returns array of length n of integers +/- 1 (randomly chosen) that represent spins
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
  #essentially from N, we get a randomized array of spin ups and downs
  #e.g., N=2 may equal [1,-1]; N=4 may equal [1,-1,-1,1] as our initial configuration
end

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


#standard function for getting the hamilitonian of 1/r^2
#Ising Model
function get_energy(s, h, J)
  E0 = 0.0
  E1 = 0.0
  E2 = 0.0
  Threads.@threads for i=1:length(s)
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


#magnization equals the summation of spins in
#a configuration.
#magnetization per spin = <M>/N, where
#N is basically length of the configuration
function get_magnetization(config::Array{Int64})
  M = sum(config)
  return M
end

function get_susceptibility(M_list, Msq_list, kT, N)
  avg_M = 0
  #avg_M = sum(M_list)/mcsteps
  avg_Msq = sum(Msq_list)/mcsteps
  X = (1/kT)*(avg_Msq - avg_M^2)/N
  return X
end

function do_MC_Step(config, kT, J, h, M_list, Msq_list)
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
  M = get_magnetization(config)
  M_sq = M^2
  #println("M squared: ", M_sq)
  #push!(M_list, abs(M))
  push!(Msq_list, M_sq)
  #we get the magnization per spin of the updated system
end


function metropolis(config_initial::AbstractArray, kT, J, h, mcsteps)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array
  mags = Vector{Float64}()
  square_mags = Vector{Float64}()
  M = get_magnetization(config)
  push!(mags, M)
  push!(square_mags, M^2)

  accepted_states = 0
  #initial state of the configuration before a MCMC update
  E = get_energy(config, h, J)
  # hamiltonian of a particular configuration
  for i=1:mcsteps
    #For MCMC an arbitrary number steps are executed for precision
    do_MC_Step(config, kT, J, h, mags, square_mags)
    #println("MC step ", i, " at ", kT, " kT.")
  end
  return mags, square_mags
end


#CONSTANTS:

#Number of spins in initialized configuration
const N = 1000

#Interaction constant
const J = 1.0

#Initialize number of montecarlo steps
const mcsteps = 10000000

config0 = initial_config(N)

#Initialize random field(s)
h0 = zeros(Float64, N)

initkT = 0.05
iter = 0.05
finalkT = 3.0

X_list1 = Vector{Float64}() # why not initialize the exact number of elements here?
for kT = initkT:iter:finalkT
  data = metropolis(config0, kT, J, h0, mcsteps)
  X = get_susceptibility(data[1], data[2], kT, length(config0))
  push!(X_list1, X) # then update in place without growing the array
  println("At ", kT, " kT.")
end
plot(initkT:iter:finalkT, X_list1, xlabel = "Temperature", ylabel = "Susceptibility", linewidth = 2.5, label = "h = 0")
