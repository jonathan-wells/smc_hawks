#!/usr/bin/env julia

###############################################################################
## Graph type and methods
###############################################################################

"Multi-purpose Graph type used to produce homology network."
type Graph
    nodes::Set
    edges::Set
    weights::Dict
end

function add_node!(net::Graph, node::ASCIIString)
    push!(net, node)
end

function add_nodes!(net::Graph, edge::Pair)
    for node in edge[1]
        @assert typeof(node) == ASCIIString
        push!(net.nodes, node)
    end
end

function add_edge!(net::Graph, edge::Pair)
    push!(net.edges, edge[1])
    net.weights[edge[1]] = edge[2]
end

function delete_node!(net::Graph, node::ASCIIString)
    delete!(net.nodes, node)
    for edge in net.edges
        if node in edge
            delete!(net.edges, edge)
            delete!(net.weights, edge)
        end
    end
end

function degree(net::Graph, node::ASCIIString)
    indegree = 0
    outdegree = 0
    for edge in net.edges
        if node == edge[1]
            outdegree += 1
        elseif node == edge[2]
            indegree += 1
        end
    end
    return (indegree + outdegree, indegree, outdegree)
end

###############################################################################
## Homology network building
###############################################################################

"Read hhresults file and return ranked dictionary of best alignments"
function read_hhrfile(hhrfile, minlength=100.0, maxeval=0.01, minprob=15.0)
    data = open(readall, hhrfile)
    getid(line) = convert(ASCIIString, split(line, ['|', '_'])[2])
    query = getid(strip(match(r"Query\s+(\S+)", data)[1]))
    data = split(split(data, "\n\n")[2], "\n")[2:end]
    hits = Dict()
    rank = 1.0
    for line in data
        hit = getid(line)
        query == hit ? continue : align = (query, hit)
        # info = [Prob, E-val, P-val, Score, SS, Cols, QryHMM, TmpHMM, Tmplen]`
        info = split(strip(line[36:end], [' ', ')']), r"\s+|\(")
        info[9] == "" ? splice!(info, 9) : nothing
        @assert length(info) == 9
        val(ind) = parse(Float64, info[ind])
        if val(6) < minlength || val(2) > maxeval || val(1) < minprob
            continue
        elseif !(align in keys(hits))
            # Add rank and length-normalised/SS-weighted "normscore"
            normscore = (val(4) + val(5))/val(6)
            hits[align] = vcat([val(i) for i in 1:6], [rank, normscore])
            rank += 1.0
        end
    end
    ranks = [hits[k][7] for k in keys(hits)]
    if length(ranks) > 1
        rmax, rmin = max(ranks...), min(ranks...)
        for k in keys(hits)
            hits[k][7] = normrank(hits[k][7], rmax, rmin)
        end
    end
    return hits
end

"Normalise ranks between 0.01 and 1"
function normrank(r, rmax, rmin)
    return 1.0/(1.0 + 99.0*(r - rmin)/(rmax - rmin))
end

"Parse .hhr files in directory for significant alignments and return as graph."
function build_raw_network(hhrdirectory)
    files = readdir(hhrdirectory)
    filter!(f->ismatch(r"\.hhr", f), files)
    net = Graph(Set(), Set(), Dict())
    for hhrfile in files
        hits = read_hhrfile(hhrdirectory*"/"*hhrfile)
        for align in hits
            add_nodes!(net, align)
            add_edge!(net, align)
        end
    end
    return net
end

"Converts graph to an undirected version"
function dedirect!(net::Graph)
    reversed = Set([reverse(edge) for edge in net.edges])
    net.edges = intersect(net.edges, reversed)
    for node in net.nodes
        degree(net, node)[1] <= 1 ? delete_node!(net, node): nothing
    end
end

"""
Return the undirected weighted edge, based on the average of ranks.
Possible attributes are:
    Rank, Prob, E-val, P-val, Score, Sec. struc. score (SS), Cols
"""
function get_mutual_attribute(net::Graph, node1, node2, attribute)
    attrdict = Dict("Rank" => 7, "Prob" => 1, "E-val" => 2, "P-val" => 3,
                    "Score" => 4, "SS" => 5, "Cols" => 6, "normscore" => 8)
    geomean(vec) = prod(vec)^(1/length(vec))
    edge = (node1, node2)
    i = attrdict[attribute]
    mutualattr = geomean([net.weights[edge][i], net.weights[reverse(edge)][i]])
    return mutualattr
end

"Returns undirected network with normalised ranks and score."
function build_final_network(net::Graph)
    dedirect!(net)
    final = Graph(net.nodes, Set(), Dict())
    visited = Set()
    for edge in net.edges
        edge in visited ? continue : union!(visited, Set([edge,reverse(edge)]))
        weights = []
        for a in ["Rank", "Prob", "E-val", "P-val", "Score", "SS", "normscore"]
            push!(weights, get_mutual_attribute(net, edge[1], edge[2], a))
        end
        add_edge!(final, Pair(edge, weights))
    end
    return final
end

function write_network(net::Graph, outfilename)
    outfile = open(outfilename, "w")
    head = "source\ttarget\trank\tprob\teval\tpval\tscore\tss\tnscore\n"
    write(outfile, head)
    for edge in net.edges
        weights = join(net.weights[edge], "\t")
        line = join([edge[1], edge[2], weights], '\t')
        write(outfile, line*"\n")
    end
    close(outfile)
end

"Given a network file with uniprot IDs as labels, converts to HGNC gene labels"
function remap_nodenames!(filename, nodefile)
    nodes = open(readlines, nodefile)[2:end]
    nodemap = Dict([Pair(split(strip(line), r"\s+")...) for line in nodes])
    infile = open(readlines, filename)
    outfile = open(filename, "w")
    write(outfile, infile[1])
    for line in infile[2:end]
        line = split(line, "\t")
        genepair = [get(nodemap, line[1], line[1]),
                    get(nodemap, line[2], line[2])]
        line = vcat(genepair, line[3:end])
        write(outfile, join(line, "\t"))
    end
    close(outfile)
end

"""
Command line arguments:
    1. hhrfolder - path to folder containing hhsearch results from all
                   proteins of interest.
    2. nodefile - path to file mapping uniprot ids to HGNC gene ids
    3. outfile - path to outfile
"""
function main()
    rnet = build_raw_network(ARGS[1])
    fnet = build_final_network(rnet)
    write_network(fnet, ARGS[3])
    remap_nodenames(ARGS[3], ARGS[2])
end

# main()
