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
    #if i != length(s)
      #E0 += J*s[i]*s[i+1]
    #end
    for j=i:length(s)
      if j != i
        E1 += J*(s[i]-s[j])^2/(i-j)^2
      end
    end
    E2 += h[i]*s[i]
  end
  E = E0-E2
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
  mags = Vector{Float64}()
  #initializing an array of updated magnizations per spin
  push!(mags, m)
  #"pushes" an updated magnization per spin into the vector variable
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
    #in lieu of the MCMC algorithm, have a probability equation
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
        #m_1 = m
        #m_2 = m*m*N
        #x = (1/kT)*(m_2-(m_1)^2)
        #push!(X,x)
        #push suscepibility of a state into the array
        push!(mags, m)
        #the updated magnization "m" is "pushed" into the mags array
      else
        config[site] = -1*config[site]
        #if false, then configuration stays the same
      end
      #if (i % (steps/10)) == 0
        #print(i/(steps/10), ", ")
      #end
    end
    return states, mags
    #at the end of the runs, we return an array of updated magnization
    #per spin (and runs), as well as the number of states (and runs)
  end


#factor of T in front of [X]

N = 100
#number of spins in initialized configuration
initkT = 0.0
#1/β
J = 1.0 #make 1 -- Jonah
#our J constant for the ising model hamiltonian
h = ones(N)*(0.01)
#do h = .01 and h = .02
#just run it twice and calculate magnetization
#then take the difference

#initialized random field
steps = 500000
#arbitary number of states used for MCMC
runs = 5
#initializes number of runs, similar to Random walk for comparison

config0 = initial_config(N)
X_list = Vector{Float64}()
#initialized suscepibility array per run

m_factor = 1.0
iter = 0.2
finalkT = 5.0

function last10(vec)
  x = 0
  if length(vec) <= 10
    print("data not long enough")
  else
    for i=length(vec)-9:length(vec)
      x += vec[i]
    end
  end
  avgx = x/10
  return avgx
end

function minlength(vecvec)
  x = Inf
  for i=1:length(vecvec)
    size = length(vecvec[i])
    x = min(x, size)
  end
  return Int(x)
end

function avgruns(vecvec)
  m = minlength(vecvec)
  v = zeros(m)
  for i=1:m
    x = 0
    for j=1:runs
      x += vecvec[j][i]
      v[i] = x/runs
    end
  end
  return v
end

h_vals = [0.5, 1.0, 2.0]

ms1 = Vector{Float64}()
ms2 = Vector{Float64}()
for kT = initkT:iter:finalkT
  avgms = Vector{Float64}()
  mags_copy1 = Vector{Vector{Float64}}()
  mags_copy2 = Vector{Vector{Float64}}()
  #mags_copy3 = Vector{Vector{Float64}}()
  for j = 1:2
    #h = j*ones(N)*m_factor
    h = ones(N)*h_vals[j+1]
    mags_list = Vector{Vector{Float64}}()
    #initialized magnization per spin vector per run
    states_list = Vector{Int64}()
    #initialized number of states per run
    for i=1:runs
      data = mcm_mag(config0, kT, J, h, steps)
      #goes through the MCMC algorithm given a number of runs
      push!(states_list, data[1])
      #seperates data[1] "number of states" into a specified array variable
      #using push! takes into account the number of runs, i.e., a new run "pushes" a new number of states
      #into this new variable
      push!(mags_list, data[2]) #seperates data[2] "mags" into a specified array variable
      #using push! takes into account the number of runs, i.e., a new run "pushes"
      #a new array of magnization into this new variable
    end
    if j==1
      mags_copy1 = mags_list
    elseif j==2
      mags_copy2 = mags_list
    #else
      #mags_copy3 = mags_list
    end
  end
  #bound = min(minlength(mags_copy1), minlength(mags_copy2))
  #println("bound: ", bound)
  #final_mags1 = map((x)->last10(x, bound), mags_copy1)
  #final_mags2 = map((x)->last10(x, bound), mags_copy2)
  mags1 = map((x)->last10(x), mags_copy1)
  mags2 = map((x)->last10(x), mags_copy2)
  #mags3 = map((x)->last10(x), mags_copy3)
  avg_mag1 = sum(mags1)/runs
  #push!(ms1, 1/(avg_mag1/h_vals[1]))
  avg_mag2 = sum(mags2)/runs
  #push!(ms1, 1/(avg_mag2/h_vals[2]))
  #avg_mag3 = sum(mags3)/runs
  #push!(ms1, 1/(avg_mag3/h_vals[3]))
  #mag1 = avgruns(mags_copy1)
  #mag2 = avgruns(mags_copy2)
  X = (avg_mag2-avg_mag1)/m_factor
  push!(X_list, X)
end

#plot(initkT:iter:finalkT, ms1, xlabel = "Temperature", ylabel = "Average Magnetization", linecolor = :blue)
#plot!(initkT:iter:finalkT, ms2, linecolor = :red)
#plot!(initkT:iter:finalkT, ms3, linecolor = :green)

plot(initkT:iter:finalkT, X_list, xlabel = "Temperature", ylabel = "Susceptibility")
