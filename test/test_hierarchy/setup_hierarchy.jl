################################################################################
# Sets up planetary-hierarchy index matrix using an initial condition file
# of the form "test.ic" (included in "test/" directory), or in a 2d array of
# the form file = [x ; "x,y,z,..."] where x,y,z are Int64
################################################################################
function hierarchy(file::String)
    # Reads in a separates the file into variables
    file = readdlm(file)
    nbody = file[1,2]
    bins = file[2,2]
    bins = replace(bins,","," ")
    bins = readdlm(IOBuffer(bins),Int64)
    level = length(bins)

    # checks that the hierarchy can be created
    @assert bins[1] <= nbody/2 && bins[1] <= 2^level

    # Creates an empty array of the proper size
    h = zeros(Float64,nbody,nbody)

    # Starts filling the array and returns it
    bottom_level(nbody,bins,level,h,1)
    h[end,1:end] = 1
    return h
end

function hierarchy(file::Array{Any,1})

    # Separates the initial array into variables
    nbody = file[1]
    bins = file[2]
    bins = replace(bins,","," ")
    bins = readdlm(IOBuffer(bins),Int64)
    level = length(bins)

    # checks that the hierarchy can be created
    @assert bins[1] <= nbody/2 && bins[1] <= 2^level

    # Creates an empty array of the proper size
    h = zeros(Float64,nbody,nbody)

    # Starts filling the array and returns it
    bottom_level(nbody,bins,level,h,1)
    h[end,1:end] = 1
    return h
end

################################################################################
# Sets up the first level of the hierchary.
################################################################################
function bottom_level(nbody::Int64,bins::Array{Int64,2},level::Int64,h::Array{Float64},iter::Int64)

    # Fills the very first level of the hierarchy and iterates related variables
    h[1,1] = 1
    h[1,2] = -1
    if bins[1] > 1
        for i=1:bins[1]-1
            if !isdefined(:j) || isa(j,Void)
                global j = i
            end
            if (j+3) <= nbody
                h[i+1,j+2] = 1
                h[i+1,j+3] = -1
            end
            j = +(j,2)
        end
    end
    clear!(:j)
    binsp = bins[1]
    row = binsp[1] + 1
    bodies = bins[1] * 2

    # checks whether the desired hierarchy is symmetric and calls appropriate
    # filling functions
    if 2*binsp == nbody
        symmetric(nbody,bodies,bins,binsp,row,h,iter+1)
    else
        n_level(nbody,bodies,bins,binsp,row,h,iter+1)
    end
end

################################################################################
# Fills subsequent levels of the hierarchy, recursively.
#
# *** Only works for hierarchies with a max binaries per level of 2 (unless
# symmetric). ***
################################################################################
function n_level(nbody::Int64,bodies::Int64,bins::Array{Int64,2},binsp::Int64,row::Int64,h::Array{Float64},iter::Int64)

# Series of checks to know which level to fill and which bins number to use
if iter <= length(bins)
    if nbody != bodies
        if bins[iter] == binsp
            bodies = bodies + bins[iter]
        elseif bins[iter] > binsp
            bodies = bodies + 2*(bins[iter]) - 1
        end
    elseif nbody == bodies
        if row == nbody-1
            h[row,1:row-1] = 1
            h[row,row:bodies] = -1
        end
    end

    # Looks at the previous and current bins and calls n_level again
    if nbody >= bodies
        if binsp == 1
            if bins[iter] == 1
                h[row,1:bodies-binsp] = 1
                h[row,bodies-binsp + 1:bodies] = -1
                row = row + binsp
                binsp = bins[iter]
                n_level(nbody,bodies,bins,binsp,row,h,iter+1)
            elseif bins[iter] == 2
                h[row,1:bodies-(2*binsp)-1] = 1
                h[row,bodies-(2*binsp)] = -1
                h[row+1,bodies-(2*binsp)+1] = 1
                h[row+1,bodies] = -1
                row = row + 2
                binsp = bins[iter]
                n_level(nbody,bodies,bins,binsp,row,h,iter+1)
            end
        elseif binsp == 2
            if bins[iter] == 1
                h[row,1:row-1] = 1
                h[row,row:bodies] = -1
                row = row + 1
                binsp = bins[iter]
                n_level(nbody,bodies,bins,binsp,row,h,iter+1)
            elseif bins[iter] == 2
                h[row,1:row-1] = 1
                h[row,row:bodies-2*(binsp-1)] = -1
                h[row+1,bodies-2*(binsp-1)+1] = 1
                h[row+1,bodies] = -1
                row = row + 2
                binsp = bins[iter]
                n_level(nbody,bodies,bins,binsp,row,h,iter+1)
            end
        end
    end

end
end

################################################################################
# Called if the hierarchy is symmetric (or if the hierarchy has all of the
# bodies on the bottom level). Fills the remaining rows.
#
# *** Only works for 4 and 8 body hierarchies ***
################################################################################
function symmetric(nbody::Int64,bodies::Int64,bins::Array{Int64,2},binsp::Int64,row::Int64,h::Array{Float64},iter::Int64)
    global j = 0
    while row < nbody-1
        h[row,1+j:2+j] = 1
        h[row,j+3:j+4] = -1
        j = j + 4
        row = row + 1
    end
    h[row,1:binsp] = 1
    h[row,binsp+1:nbody] = -1
    clear!(:j)
end
