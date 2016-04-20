#!/usr/bin/env julia

using DataStructures

function read_hhrfile(hhrfile, minlegth=50)
    data = open(readlines, hhrfile)
    queryseq = split(data[1], r"\s+")[2]
    hits = Dict()
    pat = Regex("\\s*[0-9]+\\s+((tr|sp)[A-Z0-9_|]+)\\s+.{1,9}\\s+(\\d+.\\d)\\s+([0-9.E-]+)\\s+[0-9.E-]+\\s+[0-9.]+\\s+[0-9.]+\\s+(\\d+)")
    for line in data
        hit = match(pat, line)
        if hit != nothing
            # matched cols should be greater than length of HEAT repeat
            if parse(Int, hit[5]) > minlength && ! (hit[1] in keys(hits))
                hits[hit[1]] = (parse(hit[3]), parse(hit[4]))
            end
        end
    end
    delete!(hits, queryseq)
    return queryseq, hits
end

function build_network(hhrdir, perc)
    files = readdir(hhrdir)
    filter!(f->ismatch(r"\.hhr", f), files)
    networkname = split(hhrdir, "/")[end]*"_network.txt"
    outfile = open(hhrdir*"/"*networkname, "w")
    for hhrfile in files
        query, hits = read_hhrfile(hhrdir*"/"*hhrfile)
        ranked = sort(collect(keys(hits)), by = x->hits[x][1], rev = true)
        # numhits = length(ranked)
        # if numhits == 0
        #     continue
        # elseif 0 < numhits <= 5
        #     numedges = numhits
        # else
        #     numedges = max(5, convert(Int, round(numhits*perc)))
        # end
        r = 1
        for hit in ranked[1:end]
            q = split(split(query, "|")[3], "_")[1]
            h = split(split(hit, "|")[3], "_")[1]
            line = join([q, h, r, hits[hit][1], hits[hit][2]], "\t")
            write(outfile, line*"\n")
            r += 1
        end
    end
    close(outfile)
end

function trim_network(network_file)
    data = open(readlines, network_file)
    data = [split(line, "\t") for line in data]
    query_nodes = Set([edge[1] for edge in data])
    hit_nodes = Set([edge[2] for edge in data])
    nodes = intersect(query_nodes, hit_nodes)
    outfile = open(split(network_file, ".txt")[1]*"_trimmed.txt", "w")
    header = join(["Query", "Hit", "Rank", "Probability", "E-Value\n"], "\t")
    write(outfile, header)
    for line in data
        if line[1] in nodes && line[2] in nodes
            write(outfile, join(line, "\t"))
        end
    end
    close(outfile)
end

function mutual_rank(network_file)


function main()
    # build_network("../hhresults/spombe", 0.1)
    # build_network("../hhresults/scerevisiae", 0.1)
    trim_network("../hhresults/scerevisiae/scerevisiae_network.txt")
    trim_network("../hhresults/spombe/spombe_network.txt")
end

main()
