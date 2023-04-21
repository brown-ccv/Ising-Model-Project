using Random
using .Threads
using Serialization


#takes in an integer n
#returns array of length n of integers +/- 1 (randomly chosen) that represent spins
function random_unit_config(n::Int, h::Int, seed::Int)
  Random.seed!(seed)
  return rand([-h, h], n)
  #essentially from N, we get a randomized array of spin ups and downs
  #e.g., N=2 may equal [1,-1]; N=4 may equal [1,-1,-1,1] as our initial configuration
end

#standard function for getting the hamilitonian of 1/r^2
#Ising Model
function get_energy(s::Array{Int}, J, N::Int)
  E1 = 0.0
  E2 = 0.0
  for i = 1:N-1
    for j = i+1:N
      E1 += (s[i] - s[j])^2 / (i - j)^2
    end
    #E2 += h[i]*s[i]
  end
  E = J * E1 / 2 - E2
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

function get_susceptibility(M_acc, Msq_acc, kT, N::Int, mcsteps::Int)
  avg_M = M_acc / mcsteps
  avg_Msq = Msq_acc / mcsteps
  X = (1 / kT) * (avg_Msq - avg_M^2) / N
  return X
end

function do_MC_Step(config, kT, J, M_acc, Msq_acc, N::Int)
  E = get_energy(config, J, N)
  for i = 1:N
    site = rand(1:N)
    config[site] = -1 * config[site]
    #attempt to update one site of the configuration
    E_new = get_energy(config, J, N)
    #look at the hamiltonian of the configuration
    #given this updated site
    dE = E_new - E
    #look at the difference between old and new hamiltonian
    prob = exp(-dE / (kT))
    #Metropolis acceptance ratio
    r = rand()

    if min(1, prob) > r
      #if true, configuration is "updated"
      #accepted_states += 1
      E += dE
    else
      config[site] = -1 * config[site]
      #if false, then configuration stays the same
    end
  end
  M = get_magnetization(config)
  return M
end


function metropolis(config_initial::Vector{Int64}, kT, J, mcsteps::Int, N)
  config = copy(config_initial)
  #starts out with a configuration of randomized N spin array
  #M = get_magnetization(config)
  mag_acc = 0.0
  sq_mag_acc = 0.0

  #accepted_states = 0
  #initial state of the configuration before a MCMC update
  # hamiltonian of a particular configuration
  for i = 1:mcsteps
    #For MCMC an arbitrary number steps are executed for precision
    m = do_MC_Step(config, kT, J, mag_acc, sq_mag_acc, N)
    mag_acc += abs(m)
    sq_mag_acc += m^2
  end
  return mag_acc, sq_mag_acc
end

function execute(kT, config, J, N, mcsteps)
  mag, mag_sq = metropolis(config, kT, J, mcsteps, N)
  X = get_susceptibility(mag, mag_sq, kT, N, mcsteps)
  return (X = X, mag = mag)
end

function main()
  kT = parse(Float64, first(ARGS))
  println("kT: ", kT)
  #CONSTANTS:
  #Number of spins in initialized configuration
  N = 20
  println("N: ", N)

  #Interaction constant
  J = 1.0

  #Initialize number of montecarlo steps
  mcsteps = 10_000
  println("Number of steps: ", mcsteps)

  #random seed
  seed = rand(1:1000)

  #RUN SIMULATION:
  #Initialize configuration
  config0 = random_unit_config(N, 1, seed)

  @time  data = execute(kT, config0, J, N, mcsteps)
  #println(X)
  #write X to a txt file
  path = joinpath(@__DIR__, "data", "data-$N-$kT.txt")
  #println(path)
  mkpath(dirname(path))
  serialize(path, data)
  #open(path, "w") do io
  #  write(io, X)
  #end
  
end


main()
