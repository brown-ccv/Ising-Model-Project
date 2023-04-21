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
  return rand([-1, 1], n)
end

function gaussian_rf(N)
  nd = Normal(0, 1)
  return rand(nd, N)
end

function get_energy(s::AbstractArray, h::AbstractArray, J)
  E1 = 0.0
  E2 = 0.0
  for i=1:length(s)-1
    for j=i+1:length(s)
      E1 += (s[i]-s[j])^2/(i-j)^2
    end
    E2 += h[i]*s[i]
  end
  E2 += last(h)*last(s)
  E = J*E1/2 - E2
  return E
end


function get_magnetization(config)
  m = sum(config)
  mag_per_spin = m/N
  return mag_per_spin
end

function get_avg_mag(m_acc, n_steps)
  return m_acc/n_steps
end

function do_MC_Step(config, kT, J, h)
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
  #M_acc += abs(M)
  #M_acc += M
  #we get the magnization per spin of the updated system
  return M
end


function metropolis(config_initial::AbstractArray, kT, J, h, mcsteps)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array
  M = 0
  # hamiltonian of a particular configuration
  #for i=1:mcsteps/2
    #For MCMC an arbitrary number steps are executed for precision
    #do_MC_Step(config, kT, J, h)
    #println("MC step ", i, " at ", kT, " kT.")
  #end
  M_acc = 0
  Msq_acc = 0
  for i=1:mcsteps
    M = do_MC_Step(config, kT, J, h)
    M_acc += M
    Msq_acc += M^2
  end
  return M_acc, Msq_acc
end


#CONSTANTS:

#Number of spins in initialized configuration
N = 60

#Interaction constant
J = 1.0

#Initialize number of montecarlo steps
mcsteps = 100000

#number of simulations
num_simulations = 1

config0 = initial_config(N)

#Initialize random field(s)
h0 = zeros(N)
h1 = 0.1*ones(N)

initkT = 0.05
iter = 0.05
finalkT = 3.0

num_points = Int((finalkT - initkT)/iter + 1)

# Equilibriate system with Monte Carlo,
# and calculate (absolute value of) avg magnetization

acc_M_list = zeros(num_points)
#acc_Msq_list = zeros(num_points)

for i=1:num_simulations
  M_list = Vector{Float64}()
  #Msq_list = Vector{Float64}()
  hg = gaussian_rf(N)
  for kT = initkT:iter:finalkT
    m_acc, msq_acc = metropolis(config0, kT, J, hg, mcsteps)
    avg_mag = get_avg_mag(abs(m_acc), mcsteps)
    #avgsq_mag = get_avg_mag(msq_acc, mcsteps)
    push!(M_list, avg_mag)
    #push!(Msq_list, avgsq_mag)
    #println("At ", kT, " kT.")
  end
  global acc_M_list += M_list
  #global acc_Msq_list += Msq_list
  println("Run ", i, " completed")
end

avg_M_list = acc_M_list/num_simulations
#avg_Msq_list = acc_Msq_list/num_simulations

plot(initkT:iter:finalkT, avg_M_list, xlabel = "Temperature (kT)", ylabel = "Magnetization", linewidth = 2.5, label = "M, Gaussian RF")
#plot!(initkT:iter:finalkT, avg_Msq_list, xlabel = "Temperature (kT)", ylabel = "Magnetization", linewidth = 2.5, label = "M^2, Gaussian RF")
