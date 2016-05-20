#!/usr/bin/env julia

using StatsBase

function clean_messy_fasta(filename, optional_must=nothing)
    data = open(readall, filename)
    data = split(data, "\n>")
    oldlen = length(data)
    data = Set(data)
    newlen = length(data)
    println(filename)
    println("Removed ", oldlen-newlen, " duplicate_sequences")
    count = 0
    for sequence in data
        if ismatch(r"(LOW_QUALITY)|(partial)|(hypothetical)|(uncharacterized)",
                   sequence)
            count += 1
            delete!(data, sequence)
            continue
        end
    end
    println("Removed additional ", count, " low qual/hypothetical/partial")
    println(length(data), " seqs remaining")
    return data
end

function relabel(data, delim="_")
    relabelled = []
    for sequence in data
        splitseq = split(sequence, "\n")
        label = split(splitseq[1], ">_")[1]
        println(label)
        # label = replace(label, r" ", delim)
        # sequence = join(splitseq[2:end], "\n")
        # push!(relabelled, join([label, sequence*"\n"], "\n"))
    end
    return relabelled
end

function printlabels(filename)
    speclist = []
    data = open(readall, filename)
    data = split(data, "\n>")
    for sequence in data
        splitseq = split(sequence, "\n")
        label = split(splitseq[1], ">_")[1]
        species = split(label, ['[', ']'])[2]
        push!(speclist, species)
    end
    return Set(speclist)
end

function filter_by_species(filename, outfilename)
    species_list = list_species(filename)
    push!(species_list, "Gorilla gorilla")
    push!(species_list, "Macaca mulatta")
    push!(species_list, "Aquila chrysaetos")
    data = open(readall, filename)
    data = split(data, "\n>")
    outfile = open(outfilename, "w")
    for item in data
        splitseq = split(item, "\n")
        label = strip(splitseq[1], ['>'])
        sequence = strip(join(splitseq[2:end], "\n"))
        sp = join(split(label, "_")[2:end], " ")
        if sp in species_list
            write(outfile, ">"*label*"\n"*sequence*"\n")
        end
    end
    close(outfile)
end

function get_longest_isoform(filename, outfilename)
    data = open(readall, filename)
    data = split(data, "\n>")
    speciesdict = Dict()
    for item in data
        splitseq = split(item, "\n")
        oldlabel = split(splitseq[1], "_>")[1]
        sequence = join(splitseq[2:end], "\n")
        species = split(oldlabel, ['[', ']'])[2]
        newspecies = join(split(species, "_")[1:2], "_")
        label = split(oldlabel, Regex(species))
        newlabel = ">"*join([split(filename, "_")[1], newspecies], "_")
        info = (">"*oldlabel, strip(sequence), length(sequence))
        if !(newspecies in keys(speciesdict))
            speciesdict[newspecies] = [info]
        else
            push!(speciesdict[newspecies], info)
        end
    end
    outfile = open(outfilename, "w")
    for species in keys(speciesdict)
        longest = sort(speciesdict[species], by=x->-x[3])
        fasta = join(longest[1][1:2], "\n")*"\n"
        write(outfile, fasta)
    end
end

function list_species(filename)
    data = open(readlines, filename)
    species = []
    for line in data
        if ismatch(r">", line)
            sp = join(split(strip(line), "_")[2:end], " ")
            push!(species, sp)
        end
    end
    species = countmap(species)
    to_use = []
    count = 0
    for sp in keys(species)
        if species[sp] >= 10
            println(sp, "\t", species[sp])
            count += species[sp]
            push!(to_use, sp)
        end
    end
    println(count)
    return to_use
end

function relabel_mrna(filename, outfilename)
    data = open(readall, filename)
    data = split(data, "\n>")
    outfile = open(outfilename, "w")
    for item in data
        splitseq = split(item, "\n")
        label = splitseq[1]
        species = match(r"[A-Z][a-z]+\s[a-z]+", label)
        gene = split(outfilename, "_")[1]
        label = ">"*gene*"_"replace(species.match, " ", "_")*"\n"
        sequence = join(splitseq[2:end], "\n")*"\n"
        write(outfile, label*sequence)
    end
    close(outfile)
end

function main()
    # list_species(ARGS[1])
    filter_by_species(ARGS[1], ARGS[2])
    # filter_by_species(ARGS[1])
    # get_longest_isoform(ARGS[1], ARGS[2])
    # relabel_mrna(ARGS[1], ARGS[2])
    # data = clean_messy_fasta(ARGS[1])
    # outfile = open(ARGS[2], "w")
    # for line in relabelled
    #     write(outfile, line)
    # end
    # close(outfile)
end

main()
