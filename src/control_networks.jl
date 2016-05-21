#!/usr/bin/env julia

include("hhr_network.jl")


function map_genes(nodefile)
    nodes = open(readlines, nodefile)[2:end]
    nodemap = Dict([Pair(split(strip(line), r"\s+")...) for line in nodes])
    return nodemap
end

function get_numhits(hhrdir, nodefile)
    nodemap = map_genes(nodefile)
    files = filter(x->ismatch(r".hhr", x), readdir(hhrdir))
    for f in files
        println(split(f, ".")[1])
    end
    get_hits(f) = Pair(split(f, ".")[1] => length(read_hhrfile(hhrdir*f)))
    numhits = Dict([get_hits(f) for f in files])
    numhits = Dict([nodemap[p] => Float64(numhits[p]) for p in keys(numhits)])
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
    csize = length(prots)
    clustfile = "clust"*string(myid())*".tmp"
    clusters = [Set(split(strip(c), "\t")) for c in open(readlines, clustfile)]
    for c in clusters
        if intersect(c, prots) == prots && length(c) <= csize
            return 1
        end
    end
    return 0
end

function cleanup_tmpfiles()
    for file in readdir(pwd())
        if ismatch(r".+[1-9]+\.tmp", file)
            rm(file)
        end
    end
end

function run_control(hhrdir, netfile, nodefile, iterations, prots)
    numhits = get_numhits(hhrdir, nodefile)
    network = load_network(netfile)
    results = ones(1:iterations)
    function doall(i)
        write_mcl_graph(numhits, network)
        cluster()
        return i*prots_in_cluster(prots)
    end
    results = pmap(doall, results)
    println(iterations, " trials, p-val = ", sum(results)/iterations)
    cleanup_tmpfile()
end



# "Q04002", "P40541", "Q04264", "Q06156", "Q06680"
# Args = hhrdir, netfile, nodefile, iterations, protlist
# run_control(ARGS[1], ARGS[2], ARGS[3], parse(Int, ARGS[4]), Set(ARGS[5:end]))
# get_numhits(ARGS[1], ARGS[2])
cleanup()
