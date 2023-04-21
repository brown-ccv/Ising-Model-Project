using Plots
using Random
using Serialization
using ThreadsX
using ArgParse

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
function get_energy(s, J, N::Int)
    E1 = 0.0
    for i = 1:N-1
        for j = i+1:N
            E1 += (s[i] - s[j])^2 / (i - j)^2
        end
    end
    E = J * E1 / 2
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

function do_MC_Step(config, kT, J, N::Int)
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


function metropolis(config0, kT, J, mcsteps::Int, N)
    # config = copy(config0)
    #starts out with a configuration of randomized N spin array
    #M = get_magnetization(config)
    mag_acc = 0.0
    sq_mag_acc = 0.0

    #accepted_states = 0
    #initial state of the configuration before a MCMC update
    # hamiltonian of a particular configuration
    for i = 1:mcsteps
        #For MCMC an arbitrary number steps are executed for precision
        m = do_MC_Step(config0, kT, J, N)
        mag_acc += abs(m)
        sq_mag_acc += m^2
    end
    return mag_acc, sq_mag_acc
end

function execute(initkT::Float64, iter::Float64, finalkT::Float64, config, J::Float64, N::Int, mcsteps::Int)
    function mapred(col, kT, J, mcsteps, N)
        mag, mag_sq = metropolis(col, kT, J, mcsteps, N)
        return get_susceptibility(mag, mag_sq, kT, N, mcsteps)
    end
    f(kT, col) = mapred(col, kT, J, mcsteps, N)

    ThreadsX.map(z -> f(z...), zip(initkT:iter:finalkT, eachcol(config)))
end

function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table settings begin
        "--nspins", "-n"
        help = "number of spins"
        arg_type = Int64
        default = 50

        "--mcsteps", "-m"
        help = "number of Monte Carlo steps"
        arg_type = Int64
        default = 100_000

        "--J"
        help = "coupling constant"
        arg_type = Float64
        default = 1.0

        "--initkT"
        help = "initial temperature"
        arg_type = Float64
        default = 0.05

        "--iter"
        help = "temperature increment"
        arg_type = Float64
        default = 0.05

        "--finalkT"
        help = "final temperature"
        arg_type = Float64
        default = 2.5

        "--seed"
        help = "random seed"
        arg_type = Int64
        default = 12345

        "--plt"
        help = "plot results"
        arg_type = Bool
        default = false
    end

    return parse_args(settings; as_symbols=true)
end

function main()
    # Parse command line arguments
    args = (; parse_commandline()...)
    for (k, v) in zip(keys(args), args)
        @info "$(lpad(k, 8)) =>  $v"
    end
    plotting = args.plt
    seed = args.seed
    Nspins = args.nspins
    mcsteps = args.mcsteps
    J = args.J
    initkT = args.initkT
    iter = args.iter
    finalkT = args.finalkT

    num_runs = round(Int, (finalkT - initkT) / iter + 1)
    @info "Doing $num_runs runs, from kT=$initkT to $finalkT, in steps of $iter"

    # Initialize configuration
    config00 = random_unit_config(Nspins, 1, seed)
    config0 = repeat(config00, 1, num_runs)

    X_list1 = execute(initkT, iter, finalkT, config0, J, Nspins, mcsteps)

    out = (kT=collect(initkT:iter:finalkT), X=X_list1, Nspins=Nspins, mcsteps=mcsteps, seed=seed)

    # write X_list1 to csv file
    path = joinpath(@__DIR__, "results--Nspins-$Nspins--mcsteps-$mcsteps", "result.jld")
    mkpath(dirname(path))
    serialize(path, out)

    if plotting
        p = plot(initkT:iter:finalkT, X_list1, xlabel="Temperature (kT)", ylabel="Magnetic Susceptibility", linewidth=2.5, label="h = 0")
        savefig(p, joinpath(dirname(path), "plot.png"))
    end
    @info "Results written to $(dirname(path))"
    @info "Run completed successfully"
end

@time main()
