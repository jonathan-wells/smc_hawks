#!/usr/bin/env julia

function loadfasta(filename)
    data = open(readlines, filename)
    fastadict = Dict()
    @assert ismatch(r">", data[1])
    label = ""
    for line in data
        if ismatch(r">", line)
            label = strip(line)
            fastadict[label] = ""
        else
            fastadict[label] *= strip(line)
        end
    end
    return fastadict
end

function reversefasta(fastadict::Dict)
    for label in keys(fastadict)
        fastadict[label] = reverse(fastadict[label])
    end
    return fastadict
end

function main(filename, outfilename)
    f = loadfasta(filename)
    r = reversefasta(f)
    outfile = open(outfilename, "w")
    for label in keys(r)
        write(outfile, label, "\n")
        write(outfile, r[label], "\n")
    end
    close(outfile)
end

main(ARGS[1], ARGS[2])
