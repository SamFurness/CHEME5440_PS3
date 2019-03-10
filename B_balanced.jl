#Sam Furness
#CHEME5440 HW3
#Written on 3/4/19

using LinearAlgebra
#To check the balance, a stoichiometric matrix is set up for each reaction
    #and flux of molecules across the boundary
Stoichiometric_matrix = Array{Float64}(undef,18,20)
Stoichiometric_matrix[1,:]=[-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];#aspartate
Stoichiometric_matrix[2,:]=[1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];#argininosuccinate
Stoichiometric_matrix[3,:]=[0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0];#fumarate
Stoichiometric_matrix[4,:]=[0 1 -1 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];#arginine
Stoichiometric_matrix[5,:]=[0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0];#urea
Stoichiometric_matrix[6,:]=[0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];#ornithine
Stoichiometric_matrix[7,:]=[0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];#carbomoyl phosphate
Stoichiometric_matrix[8,:]=[-1 0 0 1 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];#citruline
Stoichiometric_matrix[9,:]=[-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];#ATP
Stoichiometric_matrix[10,:]=[0 0 -1 0 -2 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0];#H2O
Stoichiometric_matrix[11,:]=[0 0 0 1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0];#Pi
Stoichiometric_matrix[12,:]=[1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0];#AMP
Stoichiometric_matrix[13,:]=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0]; #PPi
Stoichiometric_matrix[14,:]=[0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 -1 0 0 0 0]; #NO
Stoichiometric_matrix[15,:]=[0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 -1 0 0 0]; #O2
Stoichiometric_matrix[16,:]=[0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 1 0 0]; #H+
Stoichiometric_matrix[17,:]=[0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 0 1 0]; #NADPH
Stoichiometric_matrix[18,:]=[0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]; #NADP+


#A matrix for the element balance is also created. Rows are elements
    # and columns are metabolites
element_balance = Array{Float64}(undef,6,18)
element_balance[1,:]=[4 10 4 6 1 5 1 6 10 0 0 10 0 0 0 0 21 21];#carbon
element_balance[2,:]=[7 18 4 14 4 12 4 13 16 2 3 14 4 0 0 1 30 29];#hydrogen
element_balance[3,:]=[1 4 0 4 2 2 1 3 5 0 0 5 0 1 0 0 7 7];#nitrogen
element_balance[4,:]=[4 6 4 2 1 2 5 3 13 1 4 7 7 1 2 0 17 17];#oxygen
element_balance[5,:]=[0 0 0 0 0 0 1 0 3 0 1 1 2 0 0 0 3 3];#phosphorus
element_balance[6,:]=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3];#sulfur

#print out the E array to shows elemental balance in v_1-5r
E=element_balance*Stoichiometric_matrix
