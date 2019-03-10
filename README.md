# CHEME5440_PS3
Sam Furness
3/10/19


Part A: See CHEME5440PS#3.pdf file for explanation and CHEME5440HW3_C.jl for it's implementation

Part B: Run B_balanced.jl in Julia to see printed E matrix. This matrix proves the balanced nature of the reactions because the
  reaction fluxes (first 6 columns or v's) are 0s. See CHEME5440PS#3.pdf for more details on this atom matrix method.

Part C: Run CHEME5440HW3_C.jl which calls Flux.jl. The FBA solution to maximizing urea flux will be printed. The entire flux vector can 
  be viewed by printing the variable r.
