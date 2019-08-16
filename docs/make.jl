push!(LOAD_PATH,"../src/")
using Documenter, NbodyGradient

makedocs(
    sitename="NbodyGradient",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
	"Initial Conditions" => "init.md"
    ],
    format = Documenter.HTML(
        prettyurls = false
    )
)
#=deploydocs(
    repo = "github.com/ericagol/NbodyGradient.jl.git",
)=#
