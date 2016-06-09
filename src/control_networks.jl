#!/usr/bin/env julia

include("hhr_network.jl")
using StatsBase
using RCall

###############################################################################
## Permutation tests for networks.
###############################################################################

function get_numhits(hhrdir, nodefile)
    nodemap = map_genes(nodefile)
    files = filter(x->ismatch(r".hhr", x), readdir(hhrdir))
    get_hits(f) = Pair(split(f, ".")[1] => length(read_hhrfile(hhrdir*f)))
    numhits = Dict([get_hits(f) for f in files])
    numhits = Dict([get(nodemap, p, p) => Float64(numhits[p])
                   for p in keys(numhits)])
    return numhits
end

function load_network(netfile)
    net = [split(strip(l), "\t")[1:2] for l in open(readlines, netfile)[2:end]]
    return net
end

@everywhere function write_mcl_graph(numhits::Dict, network::Array)
    randrank(e) = string((rand(1:numhits[e[1]])*rand(1:numhits[e[2]])))
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
        if intersect(c, prots) == prots && length(c) == csize
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
    println("loading data")
    numhits = get_numhits(hhrdir, nodefile)
    network = load_network(netfile)
    results = ones(1:iterations)
    function doall(i)
        write_mcl_graph(numhits, network)
        cluster()
        return i*prots_in_cluster(prots)
    end
    println("clustering random graphs")
    results = pmap(doall, results)
    pval = sum(results)/iterations
    if pval == 0
        sf = strip(string(iterations), '0')
        printpval = rpad(string(pval), length(string(iterations)), "0")*sf
        println("ran $iterations trials, p < $printpval")
    else
        println("ran $iterations trials, p = $pval")
    end
    cleanup_tmpfiles()
end

###############################################################################
## Protein length vs Rank control
###############################################################################

function normalize_ranks(ranks)
    normranks::Vector{Float64} = []
    if length(ranks) == 0
        return normranks
    end
    rmin = min(ranks...)
    rmax = max(ranks...)
    for i in ranks
        push!(normranks, (ranks[i] - rmin)/(rmax - rmin))
    end
    return normranks
end

function length_correlation(hhrdir)
    files = filter(x->ismatch(r".hhr", x), readdir(hhrdir))
    evalues::Array{Float64}, template_lengths::Array{Float64} = [], []
    for hhrfile in files
        data = open(readall, hhrdir*hhrfile)
        getid(line) = convert(ASCIIString, split(line, ['|', '_'])[2])
        query = getid(strip(match(r"Query\s+(\S+)", data)[1]))
        data = split(split(data, "\n\n")[2], "\n")[2:end]
        protset = Set([query])
        for line in data
            hit = getid(line)
            hit in protset ? continue : push!(protset, hit)
            tmplen = parse(Float64, match(r"\((\d+)\)", line)[1])
            evalue = parse(Float64, split(strip(line[36:end]), r"\s+")[1])
            push!(evalues, evalue)
            push!(template_lengths, tmplen)
        end
    end
    correl = R"cor.test"
    results = correl(template_lengths, evalues, method = "k")
    println(results)
end

###############################################################################
## Command line stuff
###############################################################################

"""
Generate random networks
Command line Args:
    1. hhr directory - directory containing hhsearch results
    2. network file - file containing relevant homology network
    3. nodefile - file mapping uniprot ids to HGNC gene ids
    4. trials - number of random networks (trials) to generate
    5. cluster - contents of cluster being tested (e.g. PDS5B, ..., NCAPD2)
"""
function main1()
    hhrdir = ARGS[1]
    networkfile = ARGS[2]
    nodefile = ARGS[3]
    numtrials = parse(Int, ARGS[4])
    protlist = Set(ARGS[5:end])
    run_control(hhrdir, networkfile, nodefile, numtrials, protlist)
end

function main2()
    length_correlation(ARGS[1])
    # normalize_ranks([1, 2, 3, 4, 5, 6])
end

main2()
