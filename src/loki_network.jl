#!/usr/bin/env julia

function loki_network_intersect(loki_folder, euk_netfile)
    lokifiles = filter(f->ismatch(r".fasta", f), readdir(loki_folder))
    lokifiles = Set([split(f, ".")[1] for f in lokifiles])
    eukdata = open(readlines, euk_netfile)
    eukprots = []
    for line in eukdata[2:end]
        prots = split(line, "\t")[1:2]
        push!(eukprots, prots...)
    end
    eukprots = Set(eukprots)
    return intersect(eukprots, lokifiles)
end

function remove_fastafiles(dirname, protset)
    fastafiles = filter(f->ismatch(r".fasta", f), readdir(dirname))
    count = 0
    for f in fastafiles
        if !(split(f, ".")[1] in protset)
            println(f)
            rm(dirname*f)
            count += 1
        end
    end
    println("removed ", count, " files")
end

function main()
    protset = loki_network_intersect(ARGS[1], ARGS[2])
    remove_fastafiles(ARGS[1], protset)
end

main()
