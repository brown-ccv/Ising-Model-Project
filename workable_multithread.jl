import Pkg
Pkg.add("DataStructures")
Pkg.add("Shuffle")
Pkg.add("Plots")
Pkg.add("Distributions")
Pkg.add("Random")
using DataStructures
using Shuffle
using Plots
using Random
using Distributions

function initial_config(n::Int)
  config = zeros(Int, n)
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

function gaussian_rf(N)
  nd = Normal(0, 1)
  return rand(nd, N)
end

function unit_rf(N, h)
  field = zeros(N)
  for i=1:N
    field[i] = rand([-h, h])
  end
  return field
end

function get_energy(s::Vector{Int}, h::Vector{Float64}, J)
  E0 = 0.0
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
  E = -E1/2 - E2
  return E
end

function get_magnetization(config::Vector{Int})
  M = sum(config)
  return M
end

function get_susceptibility(M_list, Msq_list, kT, N)
  avg_M = 0
  avg_Msq = sum(Msq_list)/mcsteps
  X = (1/kT)*(avg_Msq - avg_M^2)/N
  return X
end

function do_MC_Step(config::Vector{Int}, kT, J, h::Vector{Float64}, M_list::Vector{Float64}, Msq_list::Vector{Float64})
  N = length(config)
  E = get_energy(config, h, J)
  for i=1:N
    site = rand(1:N)
    config[site] = -1*config[site]
    E_new = get_energy(config, h, J)
    dE = E_new - E
    prob = exp(-dE/(kT))
    r = rand(Float64)

    if min(1, prob) > r
      E += dE
    else
      config[site] = -1*config[site]
    end
  end
  M = get_magnetization(config)
  M_sq = M^2
  push!(Msq_list, M_sq)
end

function metropolis(config_initial::Vector{Int}, kT, J, h::Vector{Float64}, mcsteps)
  config = copy(config_initial)
  mags = Vector{Float64}()
  square_mags = Vector{Float64}()
  M = get_magnetization(config)
  push!(mags, M)
  push!(square_mags, M^2)

  accepted_states = 0
  E = get_energy(config, h, J)
  for i=1:mcsteps
    do_MC_Step(config, kT, J, h, mags, square_mags)
  end
  return mags, square_mags
end

N = 600
J = 1.0
mcsteps = 10000
config0= initial_config(N)
h0 = zeros(N)


initkT = 0.5
iter = 0.05
finalkT = 1.5
n_temps = Int(ceil((finalkT - initkT) / iter) + 1)  # Changed floor to ceil

X_list1 = Vector{Float64}(undef, n_temps)
@time begin
# Enable multithreading in Julia
Threads.@threads for temp_idx = 1:n_temps
  kT = initkT + (temp_idx - 1) * iter
  data = metropolis(config0, kT, J, h0, mcsteps)
  X = get_susceptibility(data[1], data[2], kT, length(config0))

  X_list1[temp_idx] = X

  @info "At $kT kT."
  @info "At $X X"
end
plot(initkT:iter:finalkT, X_list1, xlabel = "Temperature", ylabel = "Susceptibility", linewidth = 2.5, label = "h = 0")
end
