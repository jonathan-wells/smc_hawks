#!/usr/bin/env julia

using DataStructures
using StatsBase

###############################################################################
## Graph type and methods
###############################################################################

type Graph
    nodes::Set
    edges::Set
    weights::Dict
end

function add_node!(net::Graph, node::ASCIIString)
    push!(net, node)
end

function add_nodes!(net::Graph, edge::Tuple)
    for node in edge
        @assert typeof(node) == ASCIIString
        push!(net.nodes, node)
    end
end

function add_nodes!(net::Graph, edge::Pair)
    for node in edge[1]
        @assert typeof(node) == ASCIIString
        push!(net.nodes, node)
    end
end

function add_edge!(net::Graph, edge::Tuple)
    push!(net.edges, edge)
    net.weights[edge] = nothing
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

function delete_edge!(net::Graph, edge::Pair)
    delete!(net.edges, edge[1])
    delete!(net.weights, edge[1])
end

function delete_edge!(net::Graph, edge::Tuple)
    delete!(net.edges, edge)
    delete!(net.weights, edge)
end

function rename_node!(net::Graph, oldnode, newnode)
    push!(net.nodes, newnode)
    for edge in net.edges
        if oldnode in edge
            oldedge = [edge...]
            ind = find(i->i == oldnode, oldedge)[1]
            oldedge[ind] = newnode
            newedge = (oldedge...)
            push!(net.edges, newedge)
            net.weights[newedge] = net.weights[edge]
            delete!(net.edges, edge)
            delete!(net.weights, edge)
        end
    end
    pop!(net.nodes, oldnode)
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

function size(net::Graph)
    return (length(net.nodes), length(net.edges))
end

###############################################################################
## Graph building
###############################################################################

"Read hhresults file and return ranked dictionary of best alignments"
function read_hhrfile(hhrfile, minlength=100.0)
    data = open(readall, hhrfile)
    query = convert(ASCIIString, strip(match(r"Query\s+(\S+)", data)[1]))
    query = split(query, ['|', '_'])[2]
    data = split(split(data, "\n\n")[2], "\n")[2:end]
    hits = OrderedDict()
    rank = 1.0
    for line in data
        hit = split(line, ['|', '_'])[2]
        align = (convert(ASCIIString, query), convert(ASCIIString, hit))
        if align == reverse(align)
            continue
        end
        info = split(strip(line[36:end], [' ', ')']), r"\s+|\(")
        info[9] == "" ? splice!(info, 9) : nothing
        @assert length(info) == 9
        if !(align in keys(hits)) && parse(Int, info[6]) >= minlength
            # Add rank, adjusted after filtering for duplicates/thresholds
            hits[align] = push!([parse(Float64, i) for i in info[1:6]], rank)
            rank += 1.0
            # Add length-normalised and weighted SS score
            normscore = (hits[align][4] + hits[align][5])/hits[align][6]
            push!(hits[align], normscore)
        end
    end
    return hits
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

"Remove nodes with outdegree 0 from network.
Memory and time consumption of this function is a mess!
"
function dedirect(net::Graph)
    # edges = SharedArray(net.edges)
    for edge in deepcopy(net.edges)
        if reverse(edge) in net.edges
            continue
        end
        delete_edge!(net, edge)
    end
    for node in net.nodes
        if degree(net, node)[1] <= 1
            delete_node!(net, node)
        end
    end
    return net
end

function trim_network(net::Graph, maxeval=Inf, minprob=0, unidirectional=false)
    visited = Set()
    for edge in deepcopy(net.edges)
        edge in visited ? continue : union!(visited, Set([edge,reverse(edge)]))
        wt = net.weights[edge]
        if unidirectional
            if wt[2] > maxeval
                delete_edge!(net, edge)
            elseif wt[1] < minprob
                delete_edge!(net, edge)
            end
        else
            rwt = net.weights[reverse(edge)]
            if rwt[2] > maxeval || wt[2] > maxeval
                delete_edge!(net, edge)
                delete_edge!(net, reverse(edge))
                continue
            elseif rwt[1] < minprob || wt[1] < minprob
                delete_edge!(net, edge)
                delete_edge!(net, reverse(edge))
            end
        end
    end
    unidirectional ? mindegree = 0 : mindegree = 1
    for node in net.nodes
        if degree(net, node)[1] <= mindegree
            delete_node!(net, node)
        end
    end
    return net
end

"""
Return the undirected weighted edge, based on the average of ranks.
Possible attributes are:
    Rank
    Prob
    E-val
    P-val
    Score
    Secondary Struc Score (SS)
    Cols
"""
function get_mutual_attribute(net::Graph, node1, node2, attribute)
    attrdict = Dict("Rank" => 7, "Prob" => 1, "E-val" => 2, "P-val" => 3,
                    "Score" => 4, "SS" => 5, "Cols" => 6, "normscore" => 8)
    attribute == "Rank" ? mutualmean = mean : mutualmean = geomean
    edge = node1, node2
    i = attrdict[attribute]
    mutual_attr = mutualmean([net.weights[edge][i],
                             net.weights[reverse(edge)][i]])
    return mutual_attr
end

function map_genes(net::Graph, nodefile, anonymous=false)
    nodes = open(readlines, nodefile)[2:end]
    nodemap = Dict([node => node for node in net.nodes])
    if anonymous
        rlist = shuffle([string(i) for i in 1:2*length(nodes)])
        for line in nodes
            line = split(line)
            nodemap[line[1]] = pop!(rlist)
        end
    else
        for line in nodes
            line = split(strip(line))
            nodemap[line[1]] = line[2]
        end
    end
    return nodemap
end

function build_final_network(net::Graph, nodefile, anon=false)
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
    # Normalise ranks (between 1 and 100)
    ranks =  [w[1] for w in values(final.weights)]
    rmin, rmax = min(ranks...), max(ranks...)
    for edge in final.edges
        rank = final.weights[edge][1]
        normrank = 1.0/(1.0 + 99.0*(rank - rmin)/(rmax - rmin))
        push!(final.weights[edge], normrank)
    end
    # Convert node names from uniprot id to genes, or anonymised number
    nodemap = map_genes(net, nodefile, anon)
    for node in deepcopy(final.nodes)
        rename_node!(final, node, nodemap[node])
    end
    return final
end

function write_network(net::Graph, outfilename)
    outfile = open(outfilename, "w")
    head = "source\ttarget\trank\tprob\teval\tpval\tscore\tss\tnscore\tnrank\n"
    write(outfile, head)
    for edge in net.edges
        weights = join(net.weights[edge], "\t")
        line = join([edge[1], edge[2], weights], '\t')
        write(outfile, line*"\n")
    end
    close(outfile)
end


###############################################################################
## Control - Randomly generated graphs
###############################################################################

"Writes hhr_graph in supersimple --abc format for mcl"
function write_mcl_graph(hhrdirectory, networkfile)
    files = filter(x->ismatch(r".hhr", x), readdir(hhrdirectory))
    all_hits = [read_hhrfile(hhrdirectory*"/"*f) for f in files]
    all_hits = nothing
end




function main1()
    # ARGS are: hhrfolder, nodefile, outfile, anon_outfile
    @time rnet = build_raw_network(ARGS[1])
    @time rnet = dedirect(rnet)
    # @time rnet = trim_network(rnet, 0.01, 15.0)
    # @time fnet = build_final_network(rnet, ARGS[2])
    # write_network(fnet, ARGS[3])
    # fnetanon = build_final_network(rnet, ARGS[2], true)
    # write_network(fnetanon, ARGS[4])
end

function main2()
    @time write_mcl_graph("../data/hhresults/human_network/",
                    "../data/networks/human_network_long.txt")
end

# main1()
# main2()
