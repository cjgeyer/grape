This is an R package to simulate Markov graph models whose canonical
statistics are the number of nodes of each degree.  The point is that
we condition on not having any nodes of large degree this prevents
phase transitions or degeneracy or whatever you want to call it.
Simulation is Markov chain Monte Carlo.

To install using R package remotes
(https://cran.r-project.org/package=remotes)

    library(remotes)
    install_github("cjgeyer/grape", subdir = "package/grape")

