#!/usr/bin/env julia

function make_ssfile(a3mfile, outfilename)
    data = open(readlines, a3mfile)[1:6]
    struc = data[2]
    seq = data[6]
    outfile = open(outfilename, "w")
    for i in 1:length(seq)
        line = [string(x) for x in [i, seq[i], struc[i]]]
        line = join(line, "\t")*"\n"
        write(outfile, line)
    end
    close(outfile)
end

make_ssfile(ARGS[1], ARGS[2])