import DataStructures
import Shuffle
import Plots
#
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


#magnization equals the summation of spins in
#a configuration.
#magnetization per spin = <M>/N, where
#N is basically length of the configuration
function get_magnetization(config::AbstractArray)
  m = sum(config)
  mag_per_spin = m/length(config)
  return mag_per_spin
end


#DO I NEED MAGNETIC FIELD PARAMETER??? ASK IN MEETING -- Jonah
function mcm_mag(config_initial::AbstractArray, kT, J, h, steps)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array
  N = length(config)
  #number of spins
  m = get_magnetization(config)
  #magnetization per spin
  states = 1
  #initial state of the configuration before a MCMC update
  E = get_energy(config, h, J)
  # hamiltonian of a particular configuration
  m_1 = m
  #⟨𝜎i𝜎j⟩
  m_2 =m*m
  #⟨𝜎i⟩⟨𝜎j⟩
  x= m_2-m_1
  #susceptibility of a given state
  X = Vector{Float64}()
  #initialized susceptibility vector
  mags = Vector{Float64}()
  #initializing an array of updated magnizations per spin
  push!(mags, m)
  #"pushes" an updated magnization per spin into the vector variable
  push!(X, x)
  #"pushes" an updated suscepbility into the vector variable
  for i=1:steps
    #For MCMC an arbitrary number steps are executed for precision
    site = rand(1:N)
    #For MCMC instead of going through several configuration
    #for a given N, e.g., N=2 has 4 configurations, N=3 has 8 etc.
    #MCMC avoids that process through a randomized algorithm
    config[site] = -1*config[site]
    #attempt to update one site of the configuration
    E_new = get_energy(config, h, J)
    #look at the hamiltonian of the configuration
    #given this updated site
    dE = E_new - E
    #look at the difference between old and new hamiltonian
    prob = exp(-dE/(kT))
    #acceptance ratio equation
    r = rand(Float64)
    #as well as a randomized r from 0 to 1
      if min(1, prob) > r
       #if true, configuration is "updated"
        E = E_new
        #Hamiltonian is updated
        m = get_magnetization(config)
        #we get the magnization per spin of the updated system
        states += 1
        #clarifying that we are now in a new state
        m_1 = m
        m_2 = m*m
        x = m_2-m_1
        push!(X,x)
        #push suscepibility of a state into the array
        push!(mags, m)
        #the updated magnization "m" is "pushed" into the mags array
        
      else
        config[site] = -1*config[site]
        #if false, then configuration stays the same
      end
    end
    return states, mags, X
    #at the end of the steps, we return an array of updated magnization
    #per spin, as well as the number of states 
  end


N = 4
#number of spins in initialized configuration
kT = 5.0
#1/β
J = 1.0 #make 1 -- Jonah
#our J constant for the ising model hamiltonian
h = zeros(N)
#initialized random field
steps =99
#arbitary number of states used for MCMC
runs = 2
#initializes number of runs, similar to Random walk for comparison

config0 = initial_config(N)
mags_list = Vector{Vector{Float64}}()
#initialized magnization per spin vector per run
states_list = Vector{Int64}()
#initialized number of states per run
susceptibility_list = Vector{Vector{Float64}}()
#initialized suscepibility array per run



for i=1:runs
  data = mcm_mag(config0, kT, J, h, steps)
  #goes through the MCMC algorithm given a number of runs
  #println(data[1])
  # data[1]=number of states, data[2]=mags array
  push!(states_list, data[1])
  #seperates data[1] "number of states" into a specified array variable
  #using push! takes into account the number of runs, i.e., a new run "pushes" a new number of states
  #into this new variable
  push!(mags_list, data[2]) #seperates data[2] "mags" into a specified array variable
  #using push! takes into account the number of runs, i.e., a new run "pushes"
  #a new array of magnization into this new variable
  push!(susceptibility_list, data[3])
end
#println()
#println("states visited: ", states_list)
#println("magnetization of final state: ", last(mags_list))

plot()
#for i=1:runs
#  display(plot!(1:states_list[i], mags_list[i], label = "Run " * string(i)))
#end
println("sus", susceptibility_list[2])
println("mag", mags_list[2])

for i=1:runs
  display(plot!(1:states_list[i], susceptibility_list[i], label = "Run " * string(i)))
end
#
