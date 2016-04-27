#!/usr/bin/env julia

using DataStructures

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
        hit = strip(split(line[5:25])[1])
        hit = split(hit, ['|', '_'])[3]
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

"Parse .hhr files in directory for significant alignments and return as graph"
function build_network(hhrdirectory)
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

"Remove unidirectional edges from network"
function trim_network(net::Graph)
    for edge in net.edges
        if reverse(edge) in net.edges
            continue
        end
        delete_edge!(net, edge)
    end
    for node in net.nodes
        if degree(net, node)[1] == 0
            delete_node!(net, node)
        end
    end
    return net
end

"Return the undirected weight, based on the average of ranks"
function get_mutual_rank(net::Graph)
    undirected = Graph(Set(), Set(), Dict())
    visited = Set()
    for edge in net.edges
        if edge in visited
            continue
        end
        push!(visited, edge)
        push!(visited, reverse(edge))
        newrank = 1/mean([net.weights[edge][7], net.weights[reverse(edge)][7]])
        add_nodes!(undirected, edge)
        add_edge!(undirected, edge => (newrank, -log(newrank)))
    end
    return undirected
end

function write_network(net::Graph, outfilename)
    outfile = open(outfilename, "w")
    write(outfile, "source\ttarget\tmutualrank\tlogrank\n")
    for edge in net.edges
        weights = net.weights[edge]
        line = join([edge[1], edge[2], weights[1], weights[2]], '\t')
        write(outfile, line*"\n")
    end
    close(outfile)
    end

function main()
    ynet = build_network("../data/hhresults/scerevisiae_old")
    ynet = trim_network(ynet)
    ynet = get_mutual_rank(ynet)
    write_network(ynet, "scerevisae.txt")
    pnet = build_network("../data/hhresults/spombe_old")
    pnet = trim_network(pnet)
    pnet = get_mutual_rank(pnet)
    write_network(pnet, "spombe.txt")
end

function main2()
    ynet = build_network(ARGS[1])
    # println(size(ynet))
    ynet = trim_network(ynet)
    ynet = get_mutual_rank(ynet)
    write_network(ynet, ARGS[2])
end

main2()


