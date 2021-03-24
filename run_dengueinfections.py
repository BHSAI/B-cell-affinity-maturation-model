import math
import random
import func
import sys


#get command-line optionsccc
data_file = str(sys.argv[1])+ ".txt"
infection = str(sys.argv[2])
print( "Writing to data file " + data_file)


#immunization/infection parameters
equil_interval = 100.0
finfection_interval = 0.0 #days 
sinfection_interval = 365.0 #days

#antigen paramaters
ag_decay = 8.0 #12 for malaria
stimulation = 1.0	
ag_replication_v = 6.25
ag_replication0 = 6.2 
ag_replication1 = 6.2
ag_replication2 = 6.2
ag_replication3 = 6.2

#immune system parameters
Tcd4_initial = 1150
Tcd8_initial = 1050
Bcell_initial = 50000000	#see Smith & Perelson PNAS 1999
Bcell_carry = 5000		#see Kuppers et al. EMBO J 1993 
Bcell_aff = 10.0 
Tcell_carry = 5000
AB_aff = 2.5	
tau = 8.0			# time interval for constants (hrs) see Zhang et al Immune Lett. 1988; Liu et al Eur J. Immun. 1991

knB = round(random.triangular(0.0848, 0.2488, 0.1479), 4) 
kT4 = round(random.triangular(0.0496, 0.1223, 0.0735), 4) 
kT8 = round(random.triangular(0.0949, 0.2225, 0.1306), 4) 
print('knB, kT4, kT8 ', knB, kT4, kT8)

#define antigens
antigen_list = func.ReadAgFile( "./ep_file.txt" )

#set up antigen population
V0 = func.Antigen("V0", 0, antigen_list[0])
V1 = func.Antigen("V1", 0, antigen_list[1])
V2 = func.Antigen("V2", 0, antigen_list[2])
V3 = func.Antigen("V3", 0, antigen_list[3])

#set up populations
nB_gc = func.BCell("GC_B", 0, antigen_list )
nB_stim = func.BCell("Stimulated_B", 0, antigen_list )
nB_me = func.BCell("Memory_B", 0, antigen_list )
nB_pls = func.BCell("SL_Plasma_B", 0, antigen_list )
nB_pll = func.BCell("LL_Plasma_B", 0, antigen_list )
nAB = func.BCell("Antibody", 0, antigen_list )
print('set up population for nB')
nB = func.BCell("Naive_B", Bcell_initial, antigen_list )
print('set up population for Tcd4')
Tcd4 = func.Population("T_cd4", Tcd4_initial)
print('set up population for Tcd4_stim')
Tcd4_stim = func.Population("T_cd4_stim", 0)
print('set up population for Tcd4_memory')
Tcd4_me = func.Population("T_cd4_me", 0)
print('set up population for Tcd8')
Tcd8 = func.Population("T_cd8", Tcd8_initial)
nB_Tstim = func.BCell("Tstimulated_B", 0, antigen_list )
Tcd8_stim = func.Population("T_cd8_stim", 0)

AbClear = func.Population("AbClearance", 0)
T8Clear = func.Population("T8Clearance", 0)

GC_population = func.GroupPopulation("GC Bcell Population")
GC_population.add_population(nB_gc)
GC_population.add_population(nB_stim)
GC_population.add_population(nB_Tstim)

T4_population = func.GroupPopulation("CD4 T-cell Population")
T4_population.add_population(Tcd4)
T4_population.add_population(Tcd4_stim)

T8_population = func.GroupPopulation("CD8 T-cell Population")
T8_population.add_population(Tcd8)
T8_population.add_population(Tcd8_stim)

#set up population list
PopulationList = []
PopulationList.append(nB_gc)
PopulationList.append(nB_stim)
PopulationList.append(nB_me)
PopulationList.append(nB_pls)
PopulationList.append(nAB)
PopulationList.append(nB)
PopulationList.append(V0)
PopulationList.append(V1)
PopulationList.append(V2)
PopulationList.append(V3)
PopulationList.append(Tcd4)
PopulationList.append(Tcd4_stim)
PopulationList.append(Tcd4_me)
PopulationList.append(Tcd8)
PopulationList.append(nB_Tstim)
PopulationList.append(nB_pll)
PopulationList.append(Tcd8_stim)
PopulationList.append(AbClear)
PopulationList.append(T8Clear)


#### EQUILIBRATION 

#define a system of reactions
A0 = func.TotalReaction()

##Description: Naive B cell formation in bone marrow
eq3 = func.Formation("Naive B cell formation", float(tau/knB), nB)		#produces 250 cells every 108hrs, formerly 0.432 
A0.add_reaction(eq3)

##Description: Basal B cell decay rate
eq4b = func.Decay("Naive B cell decay", float(tau/108.0), nB)
A0.add_reaction(eq4b)

##Description: T CD4 cell formation
eq12a = func.Formation("T CD4 cell formation", float(tau/kT4), Tcd4)
A0.add_reaction(eq12a)

##Description: T CD4 cell decay rate
eq12b = func.Decay("T CD4 cell decay", float(tau/108.0), Tcd4)
A0.add_reaction(eq12b)

##Description: T CD8 cell formation
eq13a = func.Formation("T CD8 cell formation", float(tau/kT8), Tcd8)
A0.add_reaction(eq13a)

##Description: T CD8 cell decay rate
eq13b = func.Decay("T CD8 cell decay", float(tau/108.0), Tcd8)
A0.add_reaction(eq13b)

#carry out simulation
output = func.FileOutput( data_file, 0.1, PopulationList, antigen_list )
if (equil_interval > 0.0):
    output.start()

###Equilibrium (just naive B cell decay and formation) starting at day 0
total_time = float(equil_interval * 3.0 )
t = 0.0

while( t < total_time ):
	dt = A0.MC_TimeStep() 
	t += dt
	A0.MC_React()
	output.write( t )
    


#### FIRST INFECTION

ag_initial = 100 
Tstimulation1 = 1200.0  
Tstimulation2 = 2400.0  
Tstimulation3 = 8.0
Abclearance = 0.00025 
T8clearance = 0.000025

#define a system of reactions
A0 = func.TotalReaction()
print('Setting up first reaction', len(A0))

##Description: Naive B cell formation in bone marrow
eq3 = func.Formation("Naive B cell formation", float(tau/knB), nB)		#produces 250 cells every 108hrs, formerly 0.432 
A0.add_reaction(eq3)

##Description: Basal B cell decay rate
eq4b = func.Decay("Naive B cell decay", float(tau/108.0), nB)
A0.add_reaction(eq4b)

##Description: T CD4 cell formation
eq12a = func.Formation("T CD4 cell formation", float(tau/kT4), Tcd4)
A0.add_reaction(eq12a)

##Description: T CD4 cell decay rate
eq12b = func.Decay("T CD4 cell decay", float(tau/108.0), Tcd4)
A0.add_reaction(eq12b)

##Description: T CD8 cell formation
eq13a = func.Formation("T CD8 cell formation", float(tau/kT8), Tcd8)
A0.add_reaction(eq13a)

##Description: T CD8 cell decay rate
eq13b = func.Decay("T CD8 cell decay", float(tau/108.0), Tcd8)
A0.add_reaction(eq13b)

##Description: Circulating B cell antigen stimulation 
##Naive B cell stimulation is external to the GC
eq1a = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V0, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1b = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V1, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1c = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V2, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1d = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V3, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
A0.add_reaction(eq1a)
A0.add_reaction(eq1b)
A0.add_reaction(eq1c)
A0.add_reaction(eq1d)

##Description: Germinal Center B cell antigen stimulation 
##             binding affinity with antigen scales exponentially with Hamming Distance
eq2a = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V0, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity") #max stimulate rate half life 15min (0.25hr)	
eq2b = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V1, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")
eq2c = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V2, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")	
eq2d = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V3, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")
A0.add_reaction(eq2a)
A0.add_reaction(eq2b)
A0.add_reaction(eq2c)
A0.add_reaction(eq2d)

##Description: Additional stimulation of antigen-stimulated GC B cell by T_cd4 cell  
eq14 = func.TStimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 cell", float(tau/Tstimulation1), nB_stim, Tcd4, nB_Tstim, Tcd4_stim)	#cell cycle is tau
A0.add_reaction(eq14)

##Description: Additional stimulation of antigen-stimulated GC B cell by T_cd4 memory cell  
eq15 = func.TStimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 memory cell", float(tau/Tstimulation2), nB_stim, Tcd4_me, nB_Tstim, Tcd4_stim)	#cell cycle is tau
A0.add_reaction(eq15)

##Description: Stimulation of T_cd8 cell  
eq20 = func.T8Stimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 memory cell", float(tau/Tstimulation3), Tcd4_stim, Tcd8, Tcd8_stim)	#cell cycle is tau
A0.add_reaction(eq20)

##Description: Germinal Center B cell decay, governed by apoptosis. Modeled as a logistic function
##	       with a pre-defined GC population carrying capacity
eq4a = func.PopulationDecay("B Cell decay", float(tau/(tau+1.0)), float(tau/108.0), Bcell_carry, GC_population, nB_gc, nB_stim, nB_Tstim)	#base half life of 4.5 days (108hrs)
A0.add_reaction(eq4a)

eq19a = func.TPopulationDecay("Stimulated T Cell decay", float(tau/(tau+200.0)), float(tau/21600.0), Tcell_carry, T4_population, Tcd4_stim)	#base half life of ?
#eq19a = func.Decay("T CD4 cell decay", float(tau/108.0), Tcd4_stim)
A0.add_reaction(eq19a)

eq19b = func.TPopulationDecay("Stimulated T Cell decay", float(tau/(tau+200.0)), float(tau/21600.0), Tcell_carry, T8_population, Tcd8_stim)	#base half life of ?
#eq19b = func.Decay("T CD8 cell decay", float(tau/108.0), Tcd8_stim)
A0.add_reaction(eq19b)

##Description: Germinal Center B cell differentiation rate, modeling affinity-dependent T help
eq5a = func.BDifferentiation("B cell differentiation", float(tau/60.0), nB_stim, nB_gc, 1.0)	#cell cycle is tau
A0.add_reaction(eq5a)
eq5b = func.Differentiation("B cell differentiation", float(tau/8.0), nB_Tstim, nB_gc, nB_stim, nB_me, nB_pls, nB_pll, 1.0)	#cell cycle is tau
A0.add_reaction(eq5b)

##Description: TCD4 cell differentiation by mutation and into memory cell
eq16 = func.TDifferentiation("TCD4 cell differentiation", float(tau/15.0), Tcd4_stim, Tcd4, Tcd4_me)	#cell cycle is tau
A0.add_reaction(eq16)

##Description: TCD8 cell differentiation 
eq21 = func.T8Differentiation("TCD8 cell differentiation", float(tau/180.0), Tcd8_stim, Tcd8)	#cell cycle is tau
A0.add_reaction(eq21)

##Description: Antibody production by Plasma cells
eq6a = func.Production("Antibody Production", 1.0, nB_pls, nAB)
A0.add_reaction(eq6a)
eq6b = func.Production("Antibody Production", 0.1, nB_pll, nAB)
A0.add_reaction(eq6b)

##Description: Antibody decay rate
eq7 = func.Decay("Antibody Decay", float(tau/360.0), nAB)		#half life of 10 days (360hrs)
A0.add_reaction(eq7)

##Description: Intrinsic antigen decay
eq8a = func.Decay("Antigen Decay", float(tau/ag_decay), V0) 
eq8b = func.Decay("Antigen Decay", float(tau/ag_decay), V1) 
eq8c = func.Decay("Antigen Decay", float(tau/ag_decay), V2) 
eq8d = func.Decay("Antigen Decay", float(tau/ag_decay), V3) 
A0.add_reaction(eq8a)
A0.add_reaction(eq8b)
A0.add_reaction(eq8c)
A0.add_reaction(eq8d)

##Description: Antibody-dependent antigen clearance
eq8e = func.AbClearance("Antigen Clearance", Abclearance, V0, nAB, 10000.0, AB_aff, "clearance", AbClear)
eq8f = func.AbClearance("Antigen Clearance", Abclearance, V1, nAB, 10000.0, AB_aff, "clearance", AbClear)
eq8g = func.AbClearance("Antigen Clearance", Abclearance, V2, nAB, 10000.0, AB_aff, "clearance", AbClear) 
eq8h = func.AbClearance("Antigen Clearance", Abclearance, V3, nAB, 10000.0, AB_aff, "clearance", AbClear)
A0.add_reaction(eq8e)
A0.add_reaction(eq8f)
A0.add_reaction(eq8g)
A0.add_reaction(eq8h)

##Description: CD8-dependent antigen clearance
eq17a = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V0, Tcd8_stim, T8Clear) 
eq17b = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V1, Tcd8_stim, T8Clear) 
eq17c = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V2, Tcd8_stim, T8Clear) 
eq17d = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V3, Tcd8_stim, T8Clear) 
A0.add_reaction(eq17a)
A0.add_reaction(eq17b)
A0.add_reaction(eq17c)
A0.add_reaction(eq17d)

##Description: Plasma cell decay rate
eq9 = func.Decay("Plasma B cell decay", float(tau/72.0), nB_pls)		#half life of 3 days (72hrs)
A0.add_reaction(eq9)

##Description: Circulating memory cell stimulation rate
eq10a = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V0, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10b = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V1, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10c = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V2, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10d = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V3, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
A0.add_reaction(eq10a)
A0.add_reaction(eq10b)
A0.add_reaction(eq10c)
A0.add_reaction(eq10d)

##Description: viral replication 
eq11a = func.Replication("Viral replication", float(tau/ag_replication0), V0) 
eq11b = func.Replication("Viral replication", float(tau/ag_replication1), V1) 
eq11c = func.Replication("Viral replication", float(tau/ag_replication2), V2) 
eq11d = func.Replication("Viral replication", float(tau/ag_replication3), V3) 
A0.add_reaction(eq11a)
A0.add_reaction(eq11b)
A0.add_reaction(eq11c)
A0.add_reaction(eq11d)

###first infection
if (finfection_interval > 0): 
	if (infection == "monovalent"):
		V0.increase( ag_initial)   
	if (infection == "polyvalent"):
		V0.increase( ag_initial/4)
		V1.increase( ag_initial/4)
		V2.increase( ag_initial/4)
		V3.increase( ag_initial/4)
total_time = total_time + float(finfection_interval * 3.0) 
while( t <= total_time ):
	dt = A0.MC_TimeStep() 
	t += dt
	A0.MC_React()
	output.write( t )
    

#### SECOND INFECTION
    
ag_initial = 100 
Tstimulation1 = 1200.0  
Tstimulation2 = 2400.0  
Tstimulation3 = 8.0
Abclearance = 0.00025 
T8clearance = 0.000025

#define a system of reactions
A0 = func.TotalReaction()
print('Setting up second reaction', len(A0))

##Description: Naive B cell formation in bone marrow
eq3 = func.Formation("Naive B cell formation", float(tau/knB), nB)		#produces 250 cells every 108hrs, formerly 0.432 
A0.add_reaction(eq3)

##Description: Basal B cell decay rate
eq4b = func.Decay("Naive B cell decay", float(tau/108.0), nB)
A0.add_reaction(eq4b)

##Description: T CD4 cell formation
eq12a = func.Formation("T CD4 cell formation", float(tau/kT4), Tcd4)
A0.add_reaction(eq12a)

##Description: T CD4 cell decay rate
eq12b = func.Decay("T CD4 cell decay", float(tau/108.0), Tcd4)
A0.add_reaction(eq12b)

##Description: T CD8 cell formation
eq13a = func.Formation("T CD8 cell formation", float(tau/kT8), Tcd8)
A0.add_reaction(eq13a)

##Description: T CD8 cell decay rate
eq13b = func.Decay("T CD8 cell decay", float(tau/108.0), Tcd8)
A0.add_reaction(eq13b)

##Description: Circulating B cell antigen stimulation 
##Naive B cell stimulation is external to the GC
eq1a = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V0, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1b = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V1, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1c = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V2, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
eq1d = func.Stimulation("Free Naive B Cell Stimulation", float(tau/360.0), V3, nB, nB_gc, float(stimulation/6.0), Bcell_aff, "immunogenicity") 
A0.add_reaction(eq1a)
A0.add_reaction(eq1b)
A0.add_reaction(eq1c)
A0.add_reaction(eq1d)

##Description: Germinal Center B cell antigen stimulation 
##             binding affinity with antigen scales exponentially with Hamming Distance
eq2a = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V0, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity") #max stimulate rate half life 15min (0.25hr)	
eq2b = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V1, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")
eq2c = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V2, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")	
eq2d = func.Stimulation("GC B Cell Stimulation", float(tau/20.0), V3, nB_gc, nB_stim, float(stimulation/0.45), Bcell_aff, "immunogenicity")
A0.add_reaction(eq2a)
A0.add_reaction(eq2b)
A0.add_reaction(eq2c)
A0.add_reaction(eq2d)

##Description: Additional stimulation of antigen-stimulated GC B cell by T_cd4 cell  
eq14 = func.TStimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 cell", float(tau/Tstimulation1), nB_stim, Tcd4, nB_Tstim, Tcd4_stim)	#cell cycle is tau
A0.add_reaction(eq14)

##Description: Additional stimulation of antigen-stimulated GC B cell by T_cd4 memory cell  
eq15 = func.TStimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 memory cell", float(tau/Tstimulation2), nB_stim, Tcd4_me, nB_Tstim, Tcd4_stim)	#cell cycle is tau
A0.add_reaction(eq15)

##Description: Stimulation of T_cd8 cell  
eq20 = func.T8Stimulation("Stimulation of antigen-stimulated GC B cell by T_cd4 memory cell", float(tau/Tstimulation3), Tcd4_stim, Tcd8, Tcd8_stim)	#cell cycle is tau
A0.add_reaction(eq20)

##Description: Germinal Center B cell decay, governed by apoptosis. Modeled as a logistic function
##	       with a pre-defined GC population carrying capacity
eq4a = func.PopulationDecay("B Cell decay", float(tau/(tau+1.0)), float(tau/108.0), Bcell_carry, GC_population, nB_gc, nB_stim, nB_Tstim)	#base half life of 4.5 days (108hrs)
A0.add_reaction(eq4a)

eq19a = func.TPopulationDecay("Stimulated T Cell decay", float(tau/(tau+200.0)), float(tau/21600.0), Tcell_carry, T4_population, Tcd4_stim)	#base half life of ?
#eq19a = func.Decay("T CD4 cell decay", float(tau/108.0), Tcd4_stim)
A0.add_reaction(eq19a)

eq19b = func.TPopulationDecay("Stimulated T Cell decay", float(tau/(tau+200.0)), float(tau/21600.0), Tcell_carry, T8_population, Tcd8_stim)	#base half life of ?
#eq19b = func.Decay("T CD8 cell decay", float(tau/108.0), Tcd8_stim)
A0.add_reaction(eq19b)

##Description: Germinal Center B cell differentiation rate, modeling affinity-dependent T help
eq5a = func.BDifferentiation("B cell differentiation", float(tau/60.0), nB_stim, nB_gc, 1.0)	#cell cycle is tau
A0.add_reaction(eq5a)
eq5b = func.Differentiation("B cell differentiation", float(tau/8.0), nB_Tstim, nB_gc, nB_stim, nB_me, nB_pls, nB_pll, 1.0)	#cell cycle is tau
A0.add_reaction(eq5b)

##Description: TCD4 cell differentiation by mutation and into memory cell
eq16 = func.TDifferentiation("TCD4 cell differentiation", float(tau/15.0), Tcd4_stim, Tcd4, Tcd4_me)	#cell cycle is tau
A0.add_reaction(eq16)

##Description: TCD8 cell differentiation 
eq21 = func.T8Differentiation("TCD8 cell differentiation", float(tau/180.0), Tcd8_stim, Tcd8)	#cell cycle is tau
A0.add_reaction(eq21)

##Description: Antibody production by Plasma cells
eq6a = func.Production("Antibody Production", 1.0, nB_pls, nAB)
A0.add_reaction(eq6a)
eq6b = func.Production("Antibody Production", 0.1, nB_pll, nAB)
A0.add_reaction(eq6b)

##Description: Antibody decay rate
eq7 = func.Decay("Antibody Decay", float(tau/360.0), nAB)		#half life of 10 days (360hrs)
A0.add_reaction(eq7)

##Description: Intrinsic antigen decay
eq8a = func.Decay("Antigen Decay", float(tau/ag_decay), V0) 
eq8b = func.Decay("Antigen Decay", float(tau/ag_decay), V1) 
eq8c = func.Decay("Antigen Decay", float(tau/ag_decay), V2) 
eq8d = func.Decay("Antigen Decay", float(tau/ag_decay), V3) 
A0.add_reaction(eq8a)
A0.add_reaction(eq8b)
A0.add_reaction(eq8c)
A0.add_reaction(eq8d)

##Description: Antibody-dependent antigen clearance
eq8e = func.AbClearance("Antigen Clearance", Abclearance, V0, nAB, 10000.0, AB_aff, "clearance", AbClear)
eq8f = func.AbClearance("Antigen Clearance", Abclearance, V1, nAB, 10000.0, AB_aff, "clearance", AbClear)
eq8g = func.AbClearance("Antigen Clearance", Abclearance, V2, nAB, 10000.0, AB_aff, "clearance", AbClear) 
eq8h = func.AbClearance("Antigen Clearance", Abclearance, V3, nAB, 10000.0, AB_aff, "clearance", AbClear)
A0.add_reaction(eq8e)
A0.add_reaction(eq8f)
A0.add_reaction(eq8g)
A0.add_reaction(eq8h)

##Description: CD8-dependent antigen clearance
eq17a = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V0, Tcd8_stim, T8Clear) 
eq17b = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V1, Tcd8_stim, T8Clear) 
eq17c = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V2, Tcd8_stim, T8Clear) 
eq17d = func.TClearance("Antigen Clearance by CD8 T cell", T8clearance, V3, Tcd8_stim, T8Clear) 
A0.add_reaction(eq17a)
A0.add_reaction(eq17b)
A0.add_reaction(eq17c)
A0.add_reaction(eq17d)

##Description: Plasma cell decay rate
eq9 = func.Decay("Plasma B cell decay", float(tau/72.0), nB_pls)		#half life of 3 days (72hrs)
A0.add_reaction(eq9)

##Description: Circulating memory cell stimulation rate
eq10a = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V0, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10b = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V1, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10c = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V2, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
eq10d = func.MStimulation("Memory B cell Simulation", float(stimulation/80.0), V3, nB_me, nB_pls, nB_pll, float(tau/24.0), Bcell_aff, "immunogenicity") 
A0.add_reaction(eq10a)
A0.add_reaction(eq10b)
A0.add_reaction(eq10c)
A0.add_reaction(eq10d)

##Description: viral replication 
eq11a = func.Replication("Viral replication", float(tau/ag_replication0), V0) 
eq11b = func.Replication("Viral replication", float(tau/ag_replication1), V1) 
eq11c = func.Replication("Viral replication", float(tau/ag_replication2), V2) 
eq11d = func.Replication("Viral replication", float(tau/ag_replication3), V3) 
A0.add_reaction(eq11a)
A0.add_reaction(eq11b)
A0.add_reaction(eq11c)
A0.add_reaction(eq11d)

###second infection
if (sinfection_interval > 0): 
	if (infection == "monovalent"):
		PopulationList[6].increase( ag_initial)   
	if (infection == "polyvalent"):
		PopulationList[6].increase( ag_initial/4 )
		PopulationList[7].increase( ag_initial/4 )
		PopulationList[8].increase( ag_initial/4 )
		PopulationList[9].increase( ag_initial/4 )
total_time = total_time + float(sinfection_interval * 3.0) 
while( t <= total_time ):
	dt = A0.MC_TimeStep() 
	t += dt
	A0.MC_React()
	output.write( t )

    
output.finish()

