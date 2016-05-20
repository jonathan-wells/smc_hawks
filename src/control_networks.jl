#!/usr/bin/env julia

include("hhr_network.jl")

function get_numhits(hhrdir)
    files = filter(x->ismatch(r".hhr", x), readdir(hhrdir))
    get_hits(f) = Pair(split(f, ".")[1] => length(read_hhrfile(hhrdir*"/"*f)))
    numhits = Dict([get_hits(f) for f in files])
    numhits = Dict([prot => Float64(numhits[prot]) for prot in keys(numhits)])
    return numhits
end

function load_network(netfile)
    net = [split(strip(l), "\t")[1:2] for l in open(readlines, netfile)[2:end]]
    return net
end

@everywhere function write_mcl_graph(numhits, network)
    randrank(e) = string(0.5*(rand(1:numhits[e[1]])+rand(1:numhits[e[2]])))
    results = map(randrank, network)
    tmpfile = open("net"*string(myid())*".tmp", "w")
    for i in 1:length(network)
        write(tmpfile, join(network[i], "\t"), "\t", results[i], "\n")
    end
    close(tmpfile)
end

@everywhere function cluster()
    netfile = "net"*string(myid())*".tmp"
    outfile = "clust"*string(myid())*".tmp"
    run(`mcl $netfile --abc -I 2.5 -o $outfile -q x -V all`)
end

@everywhere function prots_in_cluster(prots::Set)
    clustfile = "clust"*string(myid())*".tmp"
    clusters = [Set(split(strip(c), "\t")) for c in open(readlines, clustfile)]
    for c in clusters
        if intersect(c, prots) == prots && length(c) <= length(prots)*2
            return 1
        end
    end
    return 0
end

function run_control(hhrdir, netfile, prots, iterations)
    numhits = get_numhits(hhrdir)
    network = load_network(netfile)
    results = ones(1:iterations)
    function doall(i)
        write_mcl_graph(numhits, network)
        cluster()
        return i*prots_in_cluster(prots)
    end
    results = pmap(doall, results)
    println(iterations, " trials, p-val = ", sum(results)/iterations)
    for i in procs()
        clustfile = "clust"*string(i)*".tmp"
        netfile = "net"*string(i)*".tmp"
        run(`rm $netfile $clustfile`)
    end
end

# "Q04002", "P40541", "Q04264", "Q06156", "Q06680"
@time run_control("../data/hhresults/scerevisiae_network",
            "../data/networks/scerevisiae_network_long_anon.txt",
            Set(["Q04002", "P40541", "Q04264", "Q06156", "Q06680"]), 500)
