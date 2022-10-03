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
  return config
end


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



function get_magnetization(config::AbstractArray)
  m = sum(config)
  mag_per_spin = m/length(config)
  return mag_per_spin
end


#DO I NEED MAGNETIC FIELD PARAMETER??? ASK IN MEETING
function mcm_mag(config_initial::AbstractArray, kT, J, h, steps)
  config = copy(config_initial)
  N = length(config)
  m = get_magnetization(config)
  states = 1
  E = get_energy(config, h, J)
  mags = Vector{Float64}()
  push!(mags, m)
  for i=1:steps
    site = rand(1:N)
    config[site] = -1*config[site]
    E_new = get_energy(config, h, J)
    dE = E_new - E
    prob = exp(-dE/(kT))
    r = rand(Float64)
      if min(1, prob) > r
        E = E_new
        m = get_magnetization(config)
        push!(mags, m)
        states += 1
      else
        config[site] = -1*config[site]
      end
    end
    return states, mags
  end


N = 1000
kT = 5.0
J = 1.0 #make 1
h = zeros(N)
steps = 50000
runs = 3

config0 = initial_config(N)
mags_list = Vector{Vector{Float64}}()
states_list = Vector{Int64}()

for i=1:runs
  data = mcm_mag(config0, kT, J, h, steps)
  push!(states_list, data[1])
  push!(mags_list, data[2])
end
println()
#println("states visited: ", states)
#println("magnetization of final state: ", last(mags))

plot()
for i=1:runs
  display(plot!(1:states_list[i], mags_list[i], label = "Run " * string(i)))
end
