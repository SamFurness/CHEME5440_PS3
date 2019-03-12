#Sam Furness
#CHEME5440 HW3
#Written on 3/4/19


include("Flux.jl")
#We are passing arrays of all the constraints lower and upper bounds
  #as well as the stoichiometric array in this FBA function

using LinearAlgebra

#Define the stoichiometrix matrix
Stoichiometric_matrix = Array{Float64,2}(undef,18,20)
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


#Define the bounds on the fluxes.

#Start by calculating the v upper bounds

#Define Kcat values
k_cat=Array{Float64,2}(undef,6,1)
k_cat[1]=203;#v1
k_cat[2]=34.5;#v2
k_cat[3]=249;#v3
k_cat[4]=88.1;#v4
k_cat[5]=13.7;#v5f
k_cat[6]=13.7;#v5r   #s^-1
E=0.01;#umol/gDW
theta=1;#assume enzymes maximally active

#Define metabolite concentrations from Park et al. and BRENDA
metabolite_uM=Array{Float64,2}(undef,18,1)
metabolite_uM[1]=0.0149223204061994*1000000;#aspartate
metabolite_uM[2]=0.0;#argininosucinate(not found)
metabolite_uM[3]=0.000485075251893972*1000000;#fumarate
metabolite_uM[4]=0.000255461883107884*1000000;#arginine
metabolite_uM[5]=0.0;#urea(not found)
metabolite_uM[6]=0.00448942660807892*1000000;#ornithine
metabolite_uM[7]=0.0;#carbamoyl phosphate(not found)
metabolite_uM[8]=0.0;#citruline(not found)
metabolite_uM[9]= 0.0000653952192061409*1000000#ATP
metabolite_uM[10]=0.0;#H2O (can't be found)
metabolite_uM[11]=0.0;#Pi (can't be found)
metabolite_uM[12]=0.0000423016355975951*1000000;#AMP
metabolite_uM[13]=0.0;#PPi (can't be found)
metabolite_uM[14]=0.0#NO (can't be fouond)
metabolite_uM[15]=0.0;#O2 (cant be found)
metabolite_uM[16]=0.0;#H+ (can't be found)
metabolite_uM[17]=0.0000653952192061409*1000000;#NADPH
metabolite_uM[18]=0.000502026004935134*1000000;#NADP+           #uM
#Convert uM to umol/gDW.
vol_cell=1e-12;#liters
mass_cell=2.3e-9;#g
water_fraction=0.798;#%
drymass_cell=(1-water_fraction)*mass_cell;
convertfactor(uM)=uM*(vol_cell)/drymass_cell;
metabolite=convertfactor(metabolite_uM)

#Define kM values from Park et al. and BRENDA
kM_uM=Array{Float64,2}(undef,11,1)
kM_uM[1]=0.000154276895439494*1000000;#aspartate6345 v1 Mus Musculus Park
kM_uM[2]=0.0001*1000000;#arginosuccinate4321 v2 HomoSapiens Brenda
kM_uM[3]=0.00154608643830104*1000000;#arginine3531 v3 HomoSapiens Park
kM_uM[4]=3.49694159840394;#arginine1141339 v5r Mus Musculus Park
kM_uM[5]=0.00013*1000000;#carbomoylphosphate2133 v4 HomoSapiens Brenda
kM_uM[6]=0.0016*1000000;#ornithine2133 v4 Sachromyces Cervisiae Park
kM_uM[7]=0.000056*1000000;#citruline6345 v1 HomoSapiens Brenda
kM_uM[8]=0.0;#citruline1141339 v5f (cannot be found)
kM_uM[9]=0.000392333154234199*1000000;#ATP6345   Mus Musculus Park
kM_uM[10]=0.0003*1000;#NADPH1141339 v5r MusMusculus Brenda
kM_uM[11]=0.0#NADP+ (cannot be found)   #uM
kM=convertfactor(kM_uM)

#Create flux bounds array. The default lower bounds are 0
Default_bounds_array = zeros(20,2);
#Assign the upper bounds of the v fluxes
Default_bounds_array[1,2]=k_cat[1]*E*theta*(metabolite[1]/(metabolite[1]+kM[1]))*(metabolite[9]/(metabolite[9]+kM[9]))*1;# Some terms assumed to be at saturation
Default_bounds_array[2,2]=k_cat[2]*E*theta*1;#Some terms assumed to be at saturation
Default_bounds_array[3,2]=k_cat[3]*E*theta*(metabolite[4]/(metabolite[4]+kM[3]));
Default_bounds_array[4,2]=k_cat[4]*E*theta*(metabolite[6]/(metabolite[6]+kM[6]))*1;#Some terms assumed to be at saturation
Default_bounds_array[5,2]=k_cat[5]*E*theta*1;#Some terms assumed to be at saturation
Default_bounds_array[6,2]=k_cat[6]*E*theta*(metabolite[4]/(metabolite[4]+kM[4]))*(metabolite[17]/(metabolite[17]+kM[10]));
#Assign the upper bounds for the b fluxes. Some of the metabolites (b_5-14) are reversible so the fluxes go in or out of the cell
Default_bounds_array[7,2]=(10*1000/3600);#converted to umol/gDW-s
Default_bounds_array[8,2]=(10*1000/3600);
Default_bounds_array[9,2]=(10*1000/3600);
Default_bounds_array[10,2]=(10*1000/3600);
Default_bounds_array[11,2]=(10*1000/3600);
Default_bounds_array[12,2]=(10*1000/3600);
Default_bounds_array[13,2]=(10*1000/3600);
Default_bounds_array[14,2]=(10*1000/3600);
Default_bounds_array[15,2]=(10*1000/3600);
Default_bounds_array[16,2]=(10*1000/3600);
Default_bounds_array[17,2]=(10*1000/3600);
Default_bounds_array[18,2]=(10*1000/3600);
Default_bounds_array[19,2]=(10*1000/3600);
Default_bounds_array[20,2]=(10*1000/3600);
#Set lower bounds on auxilliary metabolites as these reactions are reversible
Default_bounds_array[11,1]=-(10*1000/3600);
Default_bounds_array[12,1]=-(10*1000/3600);
Default_bounds_array[13,1]=-(10*1000/3600);
Default_bounds_array[14,1]=-(10*1000/3600);
Default_bounds_array[15,1]=-(10*1000/3600);
Default_bounds_array[16,1]=-(10*1000/3600);
Default_bounds_array[17,1]=-(10*1000/3600);
Default_bounds_array[18,1]=-(10*1000/3600);
Default_bounds_array[19,1]=-(10*1000/3600);
Default_bounds_array[20,1]=-(10*1000/3600);


#Create species bounds arrays (no bounds)
Species_bounds_array=zeros(18,2);

#Create objective function coefficient array. Only urea (b4) is maximized:
Objective_coefficient_array=zeros(20);
Objective_coefficient_array[10]=-1;

Min_flag = true;

#Use the FBA function
answer=calculate_optimal_flux_distribution(Stoichiometric_matrix,Default_bounds_array,Species_bounds_array,Objective_coefficient_array);
#Extract the flux vector
r=3600*answer[2]/1000; #mmol/gDW-hr
#Choose and print the urea flux b4
urea_flux=r[10]; #mmol/gDW-hr
println("Maximum Urea Flux in mmol/gDW*hr=",urea_flux)
