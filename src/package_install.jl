using Pkg


targets = ["Arpack","Combinatorics","Glob",
      "LaTeXStrings",
      "LatinHypercubeSampling",
      "LinearAlgebra",
      "Printf",
      "PyCall",
      "Random",
      "StatsBase",
      "Statistics",
      "SpecialFunctions",
      "TimerOutputs",
      "ThreadPools",
      "WignerSymbols"]
for tmp in targets
    Pkg.add(tmp)
end
