#!/usr/bin/env julia

using DataStructures
using StatsBase

function top_loki(hhrfile, minlength=50.0, maxeval=0.1, minprob=15.0)
    data = open(readall, hhrfile)
    query = strip(match(r"Query\s+(\S+)", data)[1])
    getid(line) = convert(ASCIIString, split(line, r"\s+")[3])
    data = split(split(data, "\n\n")[2], "\n")[2:end]
    hits = OrderedDict()  # Check me, might need to be ordered!
    rank = 1.0
    for line in data
        hit = getid(line)
        !(ismatch(r"sp\|", hit)) ? continue : nothing
        query == hit ? continue : align = (query, hit)
        # info = [Prob, E-val, P-val, Score, SS, Cols, QryHMM, TmpHMM, Tmplen]
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
    return hits
end

function loki_folder(hhrdir)
    files = readdir(hhrdir)
    for f in files
        !(ismatch(r".hhr", f)) ? continue : nothing
        hits = top_loki(hhrdir*"/"*f)
        count = 0
        for h in hits
            if count == 10
                break
            end
            println(join(h[1], "\t"), "\t", join(h[2], "\t"))
            count += 1
        end
    end
end

function tophits(filename)
    data = open(readlines, filename)
    prots = [split(line)[2] for line in data]
    cmap = countmap(prots)
    sprots = sort(collect(Set(prots)), by=x->-cmap[x])
    for i in sprots
        println(i, "\t", cmap[i])
    end
    println(sum(values(cmap)))
end

# loki_folder(ARGS[1])
tophits(ARGS[1])
