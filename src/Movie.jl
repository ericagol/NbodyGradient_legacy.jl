mutable struct Movie
   MOVIEPATH::String
   thin::Int64
   fps::Int64
   function Movie(MOVIEPATH::String;thin::Int64=1,fps::Int64=30)
       return new(MOVIEPATH,thin,fps)  
   end
end

#=
function (mov::Movie)(datafile::String;nstep::Int64)
   
    data = readdlm(datafile,comments=true)
    
    if !@isdefined(nstep)
        nstep = convert(Int64,length(data[:,1]))
    end
    prog = Progress(nstep)
    cbox = maximum([abs.(data[:,1:3:end])...,abs.(data[:,2:3:end])...]) * 1.25
    n = convert(Int64,length(data[1,:])/3)

    mp4(@animate(
    for i in 1:mov.thin:nstep 
        # Generates initial plot with proper limits, background, etc.
        scatter(data[:,1:3][i:i,1],data[:,1:3][i:i,2],marker=6,
                xaxis=("x",(-cbox,cbox)),
                yaxis=("y",(-cbox,cbox)),
                size=(600,600),framestyle=:none,bg=:black)
        # Adds planets to plot
        for j in 1:n-1
            scatter!(data[:,1+3*j:3+3*j][i:i,1],data[:,1+3*j:3+3*j][i:i,2],
                    marker=6)
        end
        # Increment progress meter
        ProgressMeter.next!(prog)
    end
    ), mov.MOVIEPATH, fps=mov.fps)
    return
end=#

function (mov::Movie)(datafile::String;nstep::Int64=0)
   
    data = readdlm(datafile,comments=true)
    
    if nstep==0
        nstep = convert(Int64,length(data[:,1]))
    end
    prog = Progress(nstep)
    cbox = maximum([abs.(data[:,1:3:end])...,abs.(data[:,2:3:end])...]) * 1.25
    n = convert(Int64,length(data[1,:])/3)

    mp4(@animate(
    for i in 1:mov.thin:nstep 
        # Generates initial plot with proper limits, background, etc.
        scatter(data[:,1:3][i:i,1],data[:,1:3][i:i,3],marker=6,
                xaxis=("x",(-cbox,cbox)),
                yaxis=("z",(-cbox,cbox)),
                size=(600,600),framestyle=:none,bg=:black)
        # Adds planets to plot
        for j in 1:n-1
            scatter!(data[:,1+3*j:3+3*j][i:i,1],data[:,1+3*j:3+3*j][i:i,3],
                    marker=6)
        end
        # Increment progress meter
        ProgressMeter.next!(prog)
    end
    ), mov.MOVIEPATH, fps=mov.fps)
    return
end
