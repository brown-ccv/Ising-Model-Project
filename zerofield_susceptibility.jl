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
  return rand([-1, 1], n)
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
  E1 = 0.0
  E2 = 0.0
  for i=1:length(s)-1
    for j=i+1:length(s)
      E1 += (s[i]-s[j])^2/(i-j)^2
    end
    #E2 += h[i]*s[i]
  end
  E = J*E1/2 - E2
  return E
end


#magnization equals the summation of spins in
#a configuration.
#magnetization per spin = M/N, where
#N is length of the configuration
function get_magnetization(config)
  M = sum(config)
  return M
end

function get_susceptibility(M_acc, Msq_acc, kT, N)
  avg_M = M_acc/mcsteps
  avg_Msq = Msq_acc/mcsteps
  X = (1/kT)*(avg_Msq - avg_M^2)/N
  return X
end

function do_MC_Step(config, kT, J, h, M_acc, Msq_acc)
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
  return M
  #M_acc += abs(M)
  #Msq_acc += M^2
  #return M_list, Msq_list
end


function metropolis(config_initial::AbstractArray, kT, J, h, mcsteps)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array
  #M = get_magnetization(config)
  mag_acc = 0
  sq_mag_acc = 0

  #accepted_states = 0
  #initial state of the configuration before a MCMC update
  E = get_energy(config, h, J)
  # hamiltonian of a particular configuration
  for i=1:mcsteps
    #For MCMC an arbitrary number steps are executed for precision
    m = do_MC_Step(config, kT, J, h, mag_acc, sq_mag_acc)
    mag_acc += abs(m)
    sq_mag_acc += m^2
  end
  return mag_acc, sq_mag_acc
end


#CONSTANTS:

#Number of spins in initialized configuration
N = 50

#Interaction constant
J = 1.0

#Initialize number of montecarlo steps
mcsteps = 1000000


#RUN SIMULATION:

#Initialize configuration
config0 = initial_config(N)

#Initialize random field(s)
h0 = zeros(Float32, N)

initkT = 0.05
iter = 0.05
finalkT = 2.5

X_list1 = Vector{Float64}()
for kT = initkT:iter:finalkT
  mag, mag_sq = metropolis(config0, kT, J, h0, mcsteps)
  X = get_susceptibility(mag, mag_sq, kT, N)
  push!(X_list1, X)
  println("At ", kT, " kT.")
end
plot(initkT:iter:finalkT, X_list1, xlabel = "Temperature (kT)", ylabel = "Magnetic Susceptibility", linewidth = 2.5, label = "h = 0")
