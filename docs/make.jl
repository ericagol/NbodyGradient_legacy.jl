using Documenter, NbodyGradient

makedocs(
    modules = [NbodyGradient],
    clean = false,
    sitename = "NbodyGradient.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
		"Initial Conditions" => "init.md"
    ],
    format = Documenter.HTML(
        prettyurls = false
    )
)

deploydocs(
    repo = "github.com/ericagol/NbodyGradient.jl.git"
)
