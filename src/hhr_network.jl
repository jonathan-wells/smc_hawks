#!/usr/bin/env julia

using DataStructures
using StatsBase

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

"Read hhresults file and return ranked dictionary of best alignments"
function read_hhrfile(hhrfile, minlength=100.0)
    data = open(readall, hhrfile)
    query = convert(ASCIIString, strip(match(r"Query\s+(\S+)", data)[1]))
    query = split(query, ['|', '_'])[3]
    data = split(split(data, "\n\n")[2], "\n")[2:end]
    hits = OrderedDict()
    rank = 1.0
    for line in data
        # hit = strip(split(line[5:25])[1])
        hit = split(line, ['|', '_'])[3]
        align = (convert(ASCIIString, query), convert(ASCIIString, hit))
        if align == reverse(align)
            continue
        end
        info = split(strip(line[36:end], [' ', ')']), r"\s+|\(")
        info[9] == "" ? splice!(info, 9) : nothing
        @assert length(info) == 9
        if !(align in keys(hits)) && parse(Int, info[6]) >= minlength
            hits[align] = push!([parse(Float64, i) for i in info[1:6]], rank)
            rank += 1.0
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

"Remove nodes with outdegree 0 from network."
function dedirect(net::Graph)
    for edge in deepcopy(net.edges)
        if reverse(edge) in net.edges
            continue
        end
        delete_edge!(net, edge)
    end
    for node in net.nodes
        if degree(net, node)[1] <= 1  # CARE NEEDED, CHECK ME
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
"""
function get_mutual_attribute(net::Graph, node1, node2, attribute)
    attrdict = Dict("Rank" => 7, "Prob" => 1, "E-val" => 2, "P-val" => 3)
    attribute == "Rank" ? mutualmean(x) = 1/mean(x) : mutualmean = geomean
    edge = node1, node2
    i = attrdict[attribute]
    mutual_attr = mutualmean([net.weights[edge][i],
                             net.weights[reverse(edge)][i]])
    return mutual_attr
end

function build_final_network(net::Graph)
    final = Graph(net.nodes, Set(), Dict())
    visited = Set()
    for edge in net.edges
        edge in visited ? continue : union!(visited, Set([edge,reverse(edge)]))
        weights = []
        for a in ["Rank", "Prob", "E-val", "P-val"]
            push!(weights, get_mutual_attribute(net, edge[1], edge[2], a))
        end
        add_edge!(final, Pair(edge, weights))
    end
    return final
end

function write_network(net::Graph, outfilename)
    outfile = open(outfilename, "w")
    write(outfile, "source\ttarget\trank\tprob\teval\tpval\n")
    for edge in net.edges
        weights = join(net.weights[edge], "\t")
        line = join([edge[1], edge[2], weights], '\t')
        write(outfile, line*"\n")
    end
    close(outfile)
end

function main()
    rnet = build_raw_network(ARGS[1])
    rnet = dedirect(rnet)
    rnet = trim_network(rnet, 0.01, 20)
    fnet = build_final_network(rnet)
    write_network(fnet, ARGS[2])
end

main()


