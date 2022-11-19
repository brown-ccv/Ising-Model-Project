import DataStructures
import Shuffle
import Plots
using DataStructures
using Shuffle
using Plots
using Random
using Distributions


#takes in an integer n
#returns array of length n of integers +/- 1 (randomly chosen) that represent spins
function initial_config(n::Int)
  config = zeros(n)
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

#standard function for getting the hamilitonian of 1/r^2
#Ising Model
function get_energy(s::AbstractArray, h::AbstractArray, J)
  E0 = 0.0
  E1 = 0.0
  E2 = 0.0
  for i=1:length(s)
    if i != length(s)
      E0 += J*s[i]*s[i+1]
    else
      E0 += J*s[i]*s[1]
    end
    #for j=i:length(s)
      #if j != i
        #E1 += J*(s[i]-s[j])^2/(i-j)^2
      #end
    #end
    #E2 += h[i]*s[i]
  end
  E = -E0
  return E
end


#magnization equals the summation of spins in
#a configuration.
#magnetization per spin = <M>/N, where
#N is basically length of the configuration
function get_magnetization(config::AbstractArray)
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
    println("MC step ", i, " at ", kT, " kT.")
  end
  return mags, square_mags
end


#CONSTANTS:

N = 100
#number of spins in initialized configuration
J = 1.0

h = zeros(N)

#initialized random field
mcsteps = 10000
#arbitary number of states used for MCMC
runs = 1
#initializes number of runs, similar to Random walk for comparison

config0 = initial_config(N)

initkT = 0.2
iter = 0.2
finalkT = 8.0

X_list = Vector{Float64}()
for kT = initkT:iter:finalkT
  data = metropolis(config0, kT, J, h, mcsteps)
  X = get_susceptibility(data[1], data[2], kT, length(config0))
  push!(X_list, X)
end

#plot(initkT:iter:finalkT, X_list, xlabel = "Temperature", ylabel = "Susceptibility")
plot(initkT:iter:finalkT, X_list, xlabel = "Temperature", ylabel = "Susceptibility", linewidth = 2.5, legend = false)
