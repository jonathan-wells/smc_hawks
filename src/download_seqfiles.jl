#!/usr/bin/env julia

using DataStructures
include("hhr_network.jl")

type HHRfile <: IO
    filename::ASCIIString
    data::Array

    function HHRfile(filename, data=nothing)
        if ismatch(r".+\.hhr", filename)
            data = open(readlines, filename)
            new(filename, data)
        else
            error("incorrect file extension")
        end
    end
end

function getallhits(hhpredresult::HHRfile)
    hitlist = []
    for line in hhpredresult.data
        hit = match(r"[0-9]{1,3}\s+.+\|+.+", line)
        if hit != nothing
            println(split(hit.match, r"\s+"))
            push!(hitlist, split(hit.match, r"\s+"))
        end
    end
    return hitlist
end

function download_fastas(hhrfile, out)
    hitlist = [split(line[2], "|")[2] for line in getallhits(hhrfile)]
    hitlist = OrderedSet(hitlist)
    url = "http://www.uniprot.org/uniprot/"
    for upr in hitlist
        println(upr)
        run(`wget -O "$out"/"$upr".fasta "$url$upr".fasta -a "$out"/wget.log`)
    end
end

function download_strict_fastas(hhrfile, out, prob=80.0, eval=0.1)
    hitlist = read_hhrfile(hhrfile, 0.0, eval, prob)
    url = "http://www.uniprot.org/uniprot/"
    for i in keys(hitlist)
        upr = i[2]
        println(upr)
        run(`wget -O "$out"/"$upr".fasta "$url$upr".fasta -a "$out"/wget.log`)
    end
end

function main()
    for file in readdir(ARGS[1])
        a = match(r"[a-zA-Z0-9]+.hhr", file)
        if a != nothing
            println(a.match)
            download_strict_fastas(a.match, ARGS[2])
        end
    end
end

main()
