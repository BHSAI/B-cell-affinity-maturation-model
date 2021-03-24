import math
import random
import itertools
import operator
import copy

ln2 = math.log(2)

#immune shape space parameters
gene_len = 20			#see Smith and Perelson PNAS 1999
gene_vocab = 4			#see Smith and Perelson PNAS 1999
max_dist = 7
min_dist = 4
germline_affinity = 6		

isotype_position = gene_len+1

#immune system parameters
differentiation_rate = 0.10	#see (1 - recycling rate) Oprea & Perelson J. Immunology 1997
mutation_rate = 0.1		#see Oprea & Perelson J. Immunology 1997, originally 0.1
isotype_rate = 0.10   #originally 0.0 
reverse_rate = 0.10
lethal_fraction = 0.0		

#########################################
##ANTIGEN BINDING AFFINITY FUNCTIONS 
#########################################

def BindingAffinity( phenotype, affinity_factor ):
	#right now binding affinity scales by ~1 with frequency of phenotype
	if (phenotype <= min_dist):
		return math.pow(affinity_factor, germline_affinity - (float(min_dist)))
	elif (phenotype <= max_dist):
		return math.pow(affinity_factor, germline_affinity - (float(phenotype)))
	else:
		return 0.0

#pre-calculated binding affinities for computational efficiency 
class BindingAffinityPrecalc: 
	def __init__( self ):
		self._BA_list = []
		self._aff = []
	
	def add_aff( self, value ):	
		BA_array = []
		for i in range(0, max_dist+1):
			BA_array.append( BindingAffinity(i, value) )
		self._aff.append( value )
		self._BA_list.append( BA_array )
		print("ADDING AFF VALUE " + str( value ))
		
	def aff_value( self, value ):
		try:
			i = self._aff.index( value )
			return i
		except ValueError:
			self.add_aff( value )
			return self._aff.index( value )
	
	def get_BA( self, phenotype, affinity_factor ):
		if (phenotype <= max_dist):
			BA_array = self._BA_list[self.aff_value( affinity_factor )]
			return BA_array[phenotype]
		else:
			return 0.0
			
BA = BindingAffinityPrecalc()
	
def BindingAffinity_Pre( phenotype, affinity_factor):
	return BA.get_BA( phenotype, affinity_factor)

#apparent size refers to the concentration of a given species (i.e. B cells)
#normalized by the binding affinity
def ApparentSize( value, phenotype, aff_factor ):
	if (float(phenotype) <= max_dist):
		n = float(value) * BindingAffinity_Pre( phenotype, aff_factor )
		return n
	else:
		return 0.0

#########################################
##BCR GENE FUNCTIONS
#########################################
def GeneRandom( ):
	sequence = ""
	for i in range( 0, gene_len ):
		sequence += str(random.randint( 1,gene_vocab ))
	
	sequence += 'M'	#add isotype
	return sequence

def GeneMutate( sequence ):
	mutation_position = random.randint( 1, gene_len )

	out_seq = ""
	for i in range( 0, gene_len ):
		if ( i == mutation_position ):
			out_seq += str(random.randint(1, gene_vocab))
		else:
			out_seq += sequence[i]
	out_seq += sequence[isotype_position-1]	#add isotype

	return out_seq

def IsotypeSwitch( sequence ):
	out_seq = ""
	for i in range(0, gene_len):
		out_seq += sequence[i]
	out_seq += 'G'	#add isotype
	return out_seq

def GenePhenotype( sequence, antigen_in ):
		ne = operator.ne
		epitope = antigen_in.epitope_all()
		dist = 999
		for i in range( 0, len( epitope ) ):
			hamming_dist = sum(map(ne, sequence, epitope[i].get_sequence()))
			if (hamming_dist < dist ):
				dist = hamming_dist
#		phenotype = min( max_dist + 1, dist )
		phenotype = dist        
		return phenotype

def GenePhenotypeEpitope( sequence, antigen_in ):
		ne = operator.ne
		epitope = antigen_in.epitope_all()
		dist = 999
		epitope_num = 999        
		for i in range( 0, len( epitope ) ):
			hamming_dist = sum(map(ne, sequence, epitope[i].get_sequence()))
			if (hamming_dist < dist ):
				dist = hamming_dist
				epitope_num = i                
#		phenotype = min( max_dist + 1, dist )
		phenotype = dist
		output = [phenotype, epitope_num]
		return output
		
def GeneEpitope( sequence, antigen_list ):
		ne = operator.ne
		dist = 999
		epitope_num = 999
		for antigen in antigen_list:
			epitope = antigen.epitope_all()
			for i in range( 0, len( epitope ) ):
				hamming_dist = sum(map(ne, sequence, epitope[i].get_sequence()))
				if (hamming_dist < dist ):
					dist = hamming_dist
					epitope_num = i
		return epitope_num
	
def GeneFromPhenotype( phenotype, antigen, ep_num ):	
	sequence = ""
	rand_pos = []
	for i in range(1, gene_len+1):
		rand_pos.append( i )
	random.shuffle(rand_pos)
	num_match = gene_len - phenotype
	epitope = antigen.epitope( ep_num ).get_sequence()

	for i in range(0,gene_len):
		if (rand_pos[i] <= num_match):
			sequence += epitope[i]
		else:
			new_val = epitope[i]
			while( new_val == epitope[i] ):
				new_val = str(random.randint(1,gene_vocab))
			sequence += new_val

	sequence += 'M'	#add isotype

	return sequence
	
#########################################
##IMMUNE SYSTEM COMPONENT TYPES
#########################################
#Name: Population
#Desc: A base class that describes the  population of a given immune system 
#      component. It is made up of subpopulations, each with their own genotype 
#      and size. There are two types of Populations, Antigen and BCell (see below).
class Population:
	def __init__( self, name, value ):
		self._name = name
		self._value = value
		self._type = ''

	def set_type ( self, type ):
		self._type = type

	def set_size ( self, value ):
		self._value = value

	def size( self ):
		return self._value

	def decrease(self, n ):
		self._value = self._value - n

	def increase(self, n ):
		self._value = self._value + n

	def name(self):
		return str(self._name)

#Name: GroupPopulation
#Desc: An aggregate population that is made up of individual Population objects
class GroupPopulation:
	def __init__( self, name ):
		self._name = name
		self._pop = []
		self._value = 0

	def add_population( self, pop ):
		self._pop.append(pop)

	def calc_population( self ):
		self._value = 0
		for i in range(0, len(self._pop) ):
			self._value += self._pop[i].size()

	def size( self ):
		self.calc_population()
		return self._value
		
	def length( self ):
		return(len(self._pop))
		
	def random_select( self ):
		self.calc_population()
		r = random.random()
		select = 0.0
		for i in range(0, len(self._pop)):
			select += float(self._pop[i].size())/float(self._value)
			if (r <= select ):
				return i
				break

	def decrease( self, n ):
		self._pop[ self.random_select() ].decrease(n)
		
	def sub_decrease( self, i, n ):
		self._pop[ i ]. decrease(n)
		
	def increase( self, n ):
		self._pop[ self.random_select() ].increase(n)
		
	def sub_increase( self, i, n ):
		self._pop[ i ].increase(n)

#Name: Epitope
#Desc: Defines a single epitope, which includes an epitope sequence (in immune shape space),
#      as well as its relative immunogenicity and clearance rates
class Epitope:
	def __init__( self, name, sequence, immunogenicity, clearance ):
		self._name = name
		self._sequence = sequence
		self._immunogenicity = immunogenicity 		
		self._clearance = clearance
	
	def mutate( self ):
		self._sequence = GeneMutate( self._sequence )
	
	def randomize( self ):
		self._sequence = GeneRandom()

	def set_sequence( self, seq ):
		self._sequence = seq
	
	def get_sequence( self ):
		return self._sequence

	def get_name ( self ):
		return self._name
		
	def immunogenicity( self ):
		return self._immunogenicity 
		
	def clearance( self ):
		return self._clearance

	def copy( self ):
		return copy.copy(self)
		
#Name: Antigen
#Desc: Defines an antigen population, which is one type of immune system component. It is 
#      defined by a name, a population size, and an AntigenType
class Antigen( Population ):
	def __init__( self, name, value, antigen ):
		Population.__init__( self, name, value )
		self._antigen = antigen
		
	def change_antigen( self, new_antigen ):
		self._antigen = new_antigen

	def return_antigen( self ):
		return self._antigen

#Name: AntigenType
#Desc: Defines an antigen type (such as a given antigen strain). This includes a name, a numerical 
#      antigenID, and a list of epitopes that make up the antigen
#note: in future version, antigenID should be set internally for ease of use
class AntigenType:
	def __init__( self, name, antigenID ):
		self._name = name
		self._epitope_array = []
		self._antigenID = antigenID
		
	def add_epitope( self, ep ):
		self._epitope_array.append( ep )
		
	def epitope( self, r ):
		return self._epitope_array[r]
	
	def epitope_all( self ):
		return self._epitope_array
	
	def epitope_num( self ):
		return len(self._epitope_array)
		
	def reset_epitopes( self ):
		self._epitope_array = []

	def ID ( self ):
		return self._antigenID
	
	def get_name ( self ):
		return self._name
		
#Name: Bcell 
#Desc: Defines a B cell population. It is indexed by genotype, and keeps track of the genotype, 
#      the size of the genotype population, the Epitope that genotype recognizes, and the phenotype
#      of that genotype with respect to all Antigens in the system.
class BCell( Population ):

	def add_genotype( self, gene, n ):
		
		self._genotype_list.append( gene )
		epitope = GeneEpitope(gene, self._antigen_list)
		self._epitope_list.append( epitope )
		self._size_list.append( n )
	
		for antigen in self._antigen_list:
			phenotype = GenePhenotype(gene, antigen)
		
			phenotype_list = self._phenotype_master[antigen.ID()]
			phenotype_list.append( phenotype )
			subpopulation = self._subpopulation_master[antigen.ID()]
			subpopulation[epitope][phenotype] += n

	def phenotype_size( self, phenotype, antigenID ):	
		pop_size = 0
		phenotype_list = self._phenotype_master[antigenID]
		for i in range(0, len(phenotype_list)):
			if (phenotype_list[i] == phenotype):
				pop_size += self._size_list[i]
		return pop_size
		
	def epitope_size( self, epitope ):
		pop_size = 0
		for i in range(0, len(self._epitope_list)):
			if (self._epitope_list[i] == epitope):
				pop_size += self._size_list[i]
		return pop_size

	def RealSizePhenotype( self, epitope, phenotype, type, antigenID):
		subpopulation = self._subpopulation_master[antigenID]
		n = subpopulation[epitope][phenotype]
		alpha = 0.0
	
		if (type == "immunogenicity"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).immunogenicity()
		elif (type == "clearance"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).clearance()
		else:
			print( "ERROR in APPARENTSIZEPHENOTYPE")
		
		real_size = n * alpha
		return real_size

	def ApparentSizePhenotype( self, epitope, phenotype, aff_factor, type, antigenID):
		subpopulation = self._subpopulation_master[antigenID]
		n = subpopulation[epitope][phenotype]
		
		alpha = 0.0
		if (type == "immunogenicity"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).immunogenicity()
		elif (type == "clearance"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).clearance()
		else:
			print("ERROR in APPARENTSIZEPHENOTYPE")
	
		app_size = ApparentSize( n, phenotype, aff_factor ) * alpha
		return app_size
		
	def ApparentSizeEpitope( self, epitope, aff_factor, type, antigenID):
		app_size = 0.0
		for i in range(0,8):
			app_size += self.ApparentSizePhenotype( epitope, i, aff_factor, type, antigenID )
		return app_size
		
	def ApparentSizeAll( self, aff_factor, type, virus ):
		antigenID = virus.return_antigen().ID()
		app_size = 0.0
		for i in range(0, self._antigen_list[0].epitope_num() ):
			app_size += self.ApparentSizeEpitope( i, aff_factor, type, antigenID)
		return app_size

	def select_random( self ):
		r = random.random()
		outcome = float(0.0)

		for i in range(0, len(self._genotype_list)):
			outcome += float(self._size_list[i])/float(Population.size( self ))
			if (r <= outcome):
				return str(self._genotype_list[i])
				break
		print("ERROR: COULD NOT FIND RANDOM select random " + str(self.name()) )
		
	def select_epitope_phenotype( self, epitope, phenotype, antigenID ):
		r = random.random()
		outcome = float(0.0)

		subpopulation = self._subpopulation_master[antigenID]
		phenotype_list = self._phenotype_master[antigenID]
		pop_size = subpopulation[epitope][phenotype]

		tot_size = 0	
		for i in range(0, len(self._genotype_list)):
			#correct = 0.0
			if (phenotype_list[i] == phenotype and self._epitope_list[i] == epitope ):
				#correct = 1.0
				tot_size += self._size_list[i]
				outcome += (float(self._size_list[i]) )/float(pop_size)
			if (r <= outcome):
				return str(self._genotype_list[i])
				break
		print("ERROR: COULD NOT FIND ###RANDOM select phenotype " + str(self._name) + ' ' + str(r) + ' ' +str(outcome) + ' ' + ' ' + str(tot_size) + ' ' + str(pop_size))
		
	def select_random_weighted( self, aff_factor, max_rate, type, agg_rate, antigen ):

		#construct phenotype probability
		r = random.random()
		
		base_phenotype = 999
		base_epitope = 999
		
		total_p = 0.0
		for i in range(0, antigen.epitope_num() ):
			for j in range(0,8):
				p_size = min(self.ApparentSizePhenotype( i, j, aff_factor, type, antigen.ID())*agg_rate , self.RealSizePhenotype(i, j, type, antigen.ID())*max_rate)
				total_p += p_size 
		
		found_it = 1
		outcome = 0.0
		for i in range(0, antigen.epitope_num() ):
			for j in range(0,8):
				p_size = min(self.ApparentSizePhenotype( i, j, aff_factor, type, antigen.ID())*agg_rate , self.RealSizePhenotype(i, j, type, antigen.ID())*max_rate)
				outcome += p_size/total_p
				if (r <= outcome  and found_it == 1):
					base_epitope = i
					base_phenotype = j
					found_it = 0
					
		if ( base_epitope == 999 or base_phenotype == 999):
			print("ERROR PHENOTYPE/EPITOPE NOT SELECTED IN REACT")
	
		subpopulation = self._subpopulation_master[antigen.ID()]	
		pop_size = subpopulation[base_epitope][base_phenotype]
		if (pop_size <= 0):
			print("ERROR: POPULATION SIZE FOR #PHENOTYPE " + str(base_phenotype) + " IS ZERO " + self.name() )
			
		return self.select_epitope_phenotype( base_epitope, base_phenotype, antigen.ID() )

	def diversity(self, threshold):
		line_count = 0
		for i in range(0, len(self._genotype_list)):
			if (self._size_list[i] >= threshold):
				line_count = line_count + 1

		return line_count
	
	def diversity2( self, threshold):
		num = threshold * float(Population.size( self ))
		sort_list = sorted( self._size_list, reverse = True )
		gene_count = 0
		pop_count = 0
		for i in range(0, len(self._genotype_list)):
			if (pop_count < num):
				pop_count += sort_list[i]
				gene_count += 1
		
		return gene_count	

	def genotype_increase( self, gene, n):
		Population.increase( self, n )
		try:
			i = self._genotype_list.index( gene )
			self._size_list[i] += n

			for antigen in self._antigen_list:
				subpopulation = self._subpopulation_master[antigen.ID()]
				phenotype_list = self._phenotype_master[antigen.ID()]
				subpopulation[self._epitope_list[i]][phenotype_list[i]] += n

		except ValueError:
			self.add_genotype( gene, n )
			
	def genotype_decrease( self, gene, n ):
		Population.decrease(self, n )
		
		try:
			i = self._genotype_list.index( gene )
			for antigen in self._antigen_list:
				phenotype_list = self._phenotype_master[antigen.ID()]
				subpopulation = self._subpopulation_master[antigen.ID()]
				subpopulation[self._epitope_list[i]][phenotype_list[i]] -= n
			
			if ((self._size_list[i] - n) <= 0):
				self._size_list.pop(i)
				self._genotype_list.pop(i)
				self._epitope_list.pop(i)
				for antigen in self._antigen_list:
					phenotype_list = self._phenotype_master[antigen.ID()]
					phenotype_list.pop(i)
			else:
				self._size_list[i] -= n
				
		except ValueError:
			print(str(self._name) + " "+ str(Population.size(self)) +str(gene)+" is not in the list!!")
		
	def increase( self, n ):
			
#		phenotype = 7
		phenotype = random.randint( 7, gene_len )
		r = random.random()
		
		n_antigen = len(self._antigen_list)
		n_epitope = len(self._antigen_list[0].epitope_all())

		r_antigen = random.random()
		r_epitope = random.random()

		count = 0.0
		antigen = 999
		for i in range(0,n_antigen):
			count = count + float(1.0/n_antigen) 
			if (r_antigen <= count and antigen == 999):
				antigen = i

		count = 0.0
		epitope = 999
		for i in range(0,n_epitope):
			count = count + float(1.0/n_epitope) 
			if (r_epitope <= count and epitope == 999):
				epitope = i
                
#		epitope = 999
#		if (r_epitope <= float(1.0/(n_antigen+1))):
#			epitope = 0
#		else:
#			epitope = 1
		
		base_gene = GeneFromPhenotype( phenotype, self._antigen_list[antigen], epitope )
		self.genotype_increase( base_gene, n )

	def decrease( self, n ):
		self.genotype_decrease( self.select_random(), 1 )

	def calc_crossreactivity_specificity( self, antigen_list ):
		output = []
		for i in range(0,60):
			output.append(0)        
		for i in range(0, len(self._genotype_list)):
			phen_epit1 = GenePhenotypeEpitope( self._genotype_list[i], self._antigen_list[0] )
			phen_epit2 = GenePhenotypeEpitope( self._genotype_list[i], self._antigen_list[1] )
			phen_epit3 = GenePhenotypeEpitope( self._genotype_list[i], self._antigen_list[2] )
			phen_epit4 = GenePhenotypeEpitope( self._genotype_list[i], self._antigen_list[3] )
			
			if (phen_epit1[0] <= max_dist and phen_epit2[0] <= max_dist and phen_epit3[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit2[1] == 0 and phen_epit3[1] == 0 and phen_epit4[1] == 0):
					output[0] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epit2[1] == 1 and phen_epit3[1] == 1 and phen_epit4[1] == 1):
					output[1] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit2[1] == 2 and phen_epit3[1] == 2 and phen_epit4[1] == 2):
					output[2] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit2[1] == 3 and phen_epit3[1] == 3 and phen_epit4[1] == 3):
					output[3] += self._size_list[i]
			elif (phen_epit1[0] <= max_dist and phen_epit2[0] <= max_dist and phen_epit3[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit2[1] == 0 and phen_epit3[1] == 0):
					output[4] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epit2[1] == 1 and phen_epit3[1] == 1):
					output[5] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit2[1] == 2 and phen_epit3[1] == 2):
					output[6] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit2[1] == 3 and phen_epit3[1] == 3):
					output[7] += self._size_list[i]
			elif (phen_epit1[0] <= max_dist and phen_epit2[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit2[1] == 0 and phen_epit4[1] == 0):
					coutput[8] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epit2[1] == 1 and phen_epit4[1] == 1):
					output[9] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit2[1] == 2 and phen_epit4[1] == 2):
					output[10] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit2[1] == 3 and phen_epit4[1] == 3):
					output[11] += self._size_list[i]
			elif (phen_epit1[0] <= max_dist and phen_epit3[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit3[1] == 0 and phen_epit4[1] == 0):
					output[12] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epi3[1] == 1 and phen_epit4[1] == 1):
					output[13] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit3[1] == 2 and phen_epit4[1] == 2):
					output[14] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit3[1] == 3 and phen_epit4[1] == 3):
					output[15] += self._size_list[i]                    
			elif (phen_epit2[0] <= max_dist and phen_epit3[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit2[1] == 0 and phen_epit3[1] == 0 and phen_epit4[1] == 0):
					output[16] += self._size_list[i]
				elif (phen_epit2[1] == 1 and phen_epi3[1] == 1 and phen_epit4[1] == 1):
					output[17] += self._size_list[i]
				elif (phen_epit2[1] == 2 and phen_epit3[1] == 2 and phen_epit4[1] == 2):
					output[18] += self._size_list[i]
				elif (phen_epit2[1] == 3 and phen_epit3[1] == 3 and phen_epit4[1] == 3):
					output[19] += self._size_list[i]                    
			elif (phen_epit1[0] <= max_dist and phen_epit2[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit2[1] == 0):
					output[20] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epi2[1] == 1):
					output[21] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit2[1] == 2):
					output[22] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit2[1] == 3):
					output[23] += self._size_list[i]
			elif (phen_epit1[0] <= max_dist and phen_epit3[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit3[1] == 0):
					output[24] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epi3[1] == 1):
					output[25] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit3[1] == 2):
					output[26] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit3[1] == 3):
					output[27] += self._size_list[i]
			elif (phen_epit1[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit1[1] == 0 and phen_epit4[1] == 0):
					output[28] += self._size_list[i]
				elif (phen_epit1[1] == 1 and phen_epi4[1] == 1):
					output[29] += self._size_list[i]
				elif (phen_epit1[1] == 2 and phen_epit4[1] == 2):
					output[30] += self._size_list[i]
				elif (phen_epit1[1] == 3 and phen_epit4[1] == 3):
					output[31] += self._size_list[i]
			elif (phen_epit2[0] <= max_dist and phen_epit3[0] <= max_dist):
				if (phen_epit2[1] == 0 and phen_epit3[1] == 0):
					output[32] += self._size_list[i]
				elif (phen_epit2[1] == 1 and phen_epi3[1] == 1):
					output[33] += self._size_list[i]
				elif (phen_epit2[1] == 2 and phen_epit3[1] == 2):
					output[34] += self._size_list[i]
				elif (phen_epit2[1] == 3 and phen_epit3[1] == 3):
					output[35] += self._size_list[i] 
			elif (phen_epit2[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit2[1] == 0 and phen_epit4[1] == 0):
					output[36] += self._size_list[i]
				elif (phen_epit2[1] == 1 and phen_epi4[1] == 1):
					output[37] += self._size_list[i]
				elif (phen_epit2[1] == 2 and phen_epit4[1] == 2):
					output[38] += self._size_list[i]
				elif (phen_epit2[1] == 3 and phen_epit4[1] == 3):
					output[39] += self._size_list[i]
			elif (phen_epit3[0] <= max_dist and phen_epit4[0] <= max_dist):
				if (phen_epit3[1] == 0 and phen_epit4[1] == 0):
					output[40] += self._size_list[i]
				elif (phen_epit3[1] == 1 and phen_epi4[1] == 1):
					output[41] += self._size_list[i]
				elif (phen_epit3[1] == 2 and phen_epit4[1] == 2):
					output[42] += self._size_list[i]
				elif (phen_epit3[1] == 3 and phen_epit4[1] == 3):
					output[43] += self._size_list[i]                   
			elif (phen_epit1[0] <= max_dist):
				if (phen_epit1[1] == 0):
					output[44] += self._size_list[i]
				elif (phen_epit1[1] == 1):
					output[45] += self._size_list[i]
				elif (phen_epit1[1] == 2):
					output[46] += self._size_list[i]
				elif (phen_epit1[1] == 3):
					output[47] += self._size_list[i]                    
			elif (phen_epit2[0] <= max_dist):
				if (phen_epit2[1] == 0):
					output[48] += self._size_list[i]
				elif (phen_epit2[1] == 1):
					output[49] += self._size_list[i]	
				elif (phen_epit2[1] == 2):
					output[50] += self._size_list[i]
				elif (phen_epit2[1] == 3):
					output[51] += self._size_list[i]
			elif (phen_epit3[0] <= max_dist):
				if (phen_epit3[1] == 0):
					output[52] += self._size_list[i]
				elif (phen_epit3[1] == 1):
					output[53] += self._size_list[i]	
				elif (phen_epit3[1] == 2):
					output[54] += self._size_list[i]
				elif (phen_epit3[1] == 3):
					output[55] += self._size_list[i]
			elif (phen_epit4[0] <= max_dist):
				if (phen_epit4[1] == 0):
					output[56] += self._size_list[i]
				elif (phen_epit4[1] == 1):
					output[57] += self._size_list[i]	
				elif (phen_epit4[1] == 2):
					output[58] += self._size_list[i]
				elif (phen_epit4[1] == 3):
					output[59] += self._size_list[i]                   
		return output

	def calc_crossreactivity( self, antigen_list ):
		output = []
		for i in range(0,16):
			output.append(0)
		for i in range(0, len(self._genotype_list)):
			phenotype1 = GenePhenotype( self._genotype_list[i], self._antigen_list[0] )
			phenotype2 = GenePhenotype( self._genotype_list[i], self._antigen_list[1] )
			phenotype3 = GenePhenotype( self._genotype_list[i], self._antigen_list[2] )
			phenotype4 = GenePhenotype( self._genotype_list[i], self._antigen_list[3] )

			if (phenotype1 <= max_dist and phenotype2 <= max_dist and phenotype3 <= max_dist and phenotype4 <= max_dist):
				output[0] += self._size_list[i]
			elif (phenotype1 <= max_dist and phenotype2 <= max_dist and phenotype3 <= max_dist):
				output[1] += self._size_list[i]
			elif (phenotype1 <= max_dist and phenotype2 <= max_dist and phenotype4 <= max_dist):
				output[2] += self._size_list[i]
			elif (phenotype1 <= max_dist and phenotype3 <= max_dist and phenotype4 <= max_dist):
				output[3] += self._size_list[i]
			elif (phenotype2 <= max_dist and phenotype3 <= max_dist and phenotype4 <= max_dist):
				output[4] += self._size_list[i]              
			elif (phenotype1 <= max_dist and phenotype2 <= max_dist):
				output[5] += self._size_list[i]
			elif (phenotype1 <= max_dist and phenotype3 <= max_dist):
				output[6] += self._size_list[i]
			elif (phenotype1 <= max_dist and phenotype4 <= max_dist):
				output[7] += self._size_list[i]
			elif (phenotype2 <= max_dist and phenotype3 <= max_dist):
				output[8] += self._size_list[i]
			elif (phenotype2 <= max_dist and phenotype4 <= max_dist):
				output[9] += self._size_list[i]
			elif (phenotype3 <= max_dist and phenotype4 <= max_dist):
				output[10] += self._size_list[i]                
			elif (phenotype1 <= max_dist):
				output[11] += self._size_list[i]
			elif (phenotype2 <= max_dist):
				output[12] +=self._size_list[i]
			elif (phenotype3 <= max_dist):
				output[13] +=self._size_list[i]
			elif (phenotype4 <= max_dist):
				output[14] +=self._size_list[i]                
			else:
				output[15] +=self._size_list[i]
		return output

	def calc_transcend( self, antigen_list ):
		output = []
		for i in range(0,len(antigen_list)+1):
			output.append(0)
		
		for i in range(0, len(self._genotype_list)):
			strain_num = 0
			#for antigen in self._antigen_list:
			for j in range(0, len(antigen_list)-1):
				phen = GenePhenotype( self._genotype_list[i], self._antigen_list[j] )
				if (phen <= max_dist):
					strain_num += 1
			output[strain_num] += self._size_list[i]
		return output

	#we need to work out depletion and neutralization 
	def calc_neutralization( self, antigen_list ):
		spec1_num = 0
		spec2_num = 0
		spec3_num = 0
		spec4_num = 0        
		for i in range(0, len(self._genotype_list)):
			phenotype1 = GenePhenotype( self._genotype_list[i], self._antigen_list[0] )
			phenotype2 = GenePhenotype( self._genotype_list[i], self._antigen_list[1] )
			phenotype3 = GenePhenotype( self._genotype_list[i], self._antigen_list[2] )
			phenotype4 = GenePhenotype( self._genotype_list[i], self._antigen_list[3] )
			
			if (phenotype1 <= max_dist):
				spec1_num += self._size_list[i]*BindingAffinity_Pre(phenotype1, 2.5)
			if (phenotype2 <= max_dist):
				spec2_num += self._size_list[i]*BindingAffinity_Pre(phenotype2, 2.5)
			if (phenotype3 <= max_dist):
				spec3_num += self._size_list[i]*BindingAffinity_Pre(phenotype3, 2.5)
			if (phenotype4 <= max_dist):
				spec4_num += self._size_list[i]*BindingAffinity_Pre(phenotype4, 2.5)
                
		output = [spec1_num, spec2_num, spec3_num, spec4_num]
		return output

	def calc_cross(self, antigen1, epitope):
		subpopulation = self._subpopulation_master[antigen1.ID()]
		pop = [0,0,0,0]
		for phenotype in range(0, 5):
			pop[0] = pop[0] + subpopulation[epitope][phenotype]
		pop[1] = subpopulation[epitope][5]
		pop[2] = subpopulation[epitope][6]
		pop[3] = subpopulation[epitope][7]
		return pop

	def calc_isotype( self ):
		g_num = 0
		for i in range(0, len(self._genotype_list)):
			gene = self._genotype_list[i]
			if (gene[isotype_position - 1] == 'G'):
				g_num += self._size_list[i]
		m_num = self.size() - g_num
		output = [m_num, g_num]
		return output
        
	#function for randomly generating the starting population of naive B cells
	def generate_population( self, value):
		for antigen in self._antigen_list:
			n_epitope = len(antigen.epitope_all())
			for j in range(0, n_epitope):
				for phenotype in range(7,8):
					current_pop = 0
					pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))/n_epitope)
					while (current_pop < pop):
						current_pop += 1
						base_gene = GeneFromPhenotype( phenotype, antigen, j )
						min_phen = 7
						for k in range(0,len(self._antigen_list)):
							phenotype1 = GenePhenotype( base_gene, self._antigen_list[k] )
							if (phenotype1 < min_phen ):
								min_phen = phenotype1 
						if (min_phen == 7):
							self.genotype_increase(base_gene, 1)
#						print( base_gene )

		high_pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))) * len(self._antigen_list) * 2
		for i in range(0, high_pop):
			for antigen in self._antigen_list:
				n_epitope = len(antigen.epitope_all())
				for j in range(0, n_epitope):
					for phenotype in range(7,8):
						pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))/n_epitope)
						subpopulation = self._subpopulation_master[ antigen.ID() ]
						current_pop = subpopulation[j][phenotype]
						if ( subpopulation[j][phenotype] > pop ):
							del_gene = self.select_epitope_phenotype(j, phenotype, antigen.ID())
							self.genotype_decrease( del_gene, 1)

	def generate_population_new( self, value):
		for antigen in self._antigen_list:
			n_epitope = len(antigen.epitope_all())
			for j in range(0, n_epitope):
				for phenotype in range(7,20):
					current_pop = 0
					pop = int((math.floor(math.pow(10,-1*(5+max_dist-7.0)) * float(value))/n_epitope)/1.0)
					while (current_pop < pop):
						current_pop += 1
						base_gene = GeneFromPhenotype( phenotype, antigen, j )
						min_phen = phenotype
						for k in range(0,len(self._antigen_list)):
							phenotype1 = GenePhenotype( base_gene, self._antigen_list[k] )
							if (phenotype1 < min_phen ):
								min_phen = phenotype1 
						if (min_phen == phenotype):
							self.genotype_increase(base_gene, 1)
#							print(current_pop, base_gene )                   
		high_pop = int(math.floor(math.pow(10,-1*(5+max_dist-7.0)) * float(value))) * len(self._antigen_list) * 2
		for i in range(0, high_pop):        
			for antigen in self._antigen_list:
				n_epitope = len(antigen.epitope_all())
				for j in range(0, n_epitope):
					for phenotype in range(7,20):
						pop = int((math.floor(math.pow(10,-1*(5+max_dist-7.0)) * float(value))/n_epitope)/1.0)
						subpopulation = self._subpopulation_master[ antigen.ID() ]
						current_pop = subpopulation[j][phenotype]
						if ( subpopulation[j][phenotype] > pop ):
							del_gene = self.select_epitope_phenotype(j, phenotype, antigen.ID())
							self.genotype_decrease( del_gene, 1)
                            
	def return_gene_list( self ):
		return self._genotype_list

	def return_size_list( self ):
		return self._size_list

	def __init__( self, name, value, antigen_list ):
		Population.__init__( self, name, 0 )
		Population.set_type(self, "bcell")
		self._genotype_list = []
		self._epitope_list = []
		self._size_list = []
		self._antigen_list = antigen_list
		epitope_num = antigen_list[0].epitope_num()
		
		#for each antigen
		self._phenotype_master = []
		self._subpopulation_master = []
		for antigen in antigen_list:	
			phenotype_list = []
			subpopulation = [[0]*20 for x in range(epitope_num)]
			
			self._phenotype_master.append(phenotype_list)	
			self._subpopulation_master.append(subpopulation)
		
		self.generate_population_new( value )
            
#########################################
## REACTIONS TYPES
#########################################
#Name: OrderZero
#Desc: Defines the base class for zero-order reactions 
class OrderZero:
    def __init__( self, name, k ):
        self._name = name
        self._k = ln2 * float( k )
        self.__rate = 0
        
    def rate( self ):
        self.__rate = self._k
        return self.__rate

#Name: OrderOne
#Desc: Defines the base class for first-order reactions
class OrderOne:
    def __init__( self, name, k, A ):
        self._name = name
        self._k = ln2*float(k)
        self._A = A

    def rate( self ):
        self.__rate = self._k * self._A.size()
        return self.__rate

#Name: OrderTwo
#Desc: Defines the base class for second-order reactions
class OrderTwo:
	def __init__( self, name, k, A, B):
		self._name = name
		self._k = ln2*float(k)
		self._A = A
		self._B = B

	def rate( self ):
		self.__rate = self._k * self._A.size() * self._B.size()
		return self.__rate

#Name: OrderTwoPhenotype
#Desc: Defines the base class for second-order reactions where the reaction rate is the 
#      result of a heterogenious population, with varying phenotypes (such as B cell stimulation)
class OrderTwoPhenotype:  
	def __init__( self, name, k, A, B, max_rate, aff_factor, type ):
		self._name = name
		self._k = ln2*float(k)
		self._A = A
		self._B = B
		self._max_rate = ln2*max_rate
		self._aff_factor = aff_factor
		BA.aff_value( aff_factor )
		self._type = type

	def rate( self ):
		self.__rate = 0.0
		if (self._B.size() > 0):
			Bcell_app = self._B.ApparentSizeAll( self._aff_factor, self._type, self._A )
			self.__rate = min(self._A.size() * Bcell_app * self._k, self._B.size() * self._max_rate)
		return self.__rate

#########################################
## IMMUNE SYSTEM REACTIONS 
#########################################
#Name: Stimulation 
#Desc: Defines the reaction for B cell stimulation
#		i.e. 	B + Ag -> B* + Ag 
class Stimulation( OrderTwoPhenotype ):
	def __init__( self, name, k, A, B, C, max_rate, aff_factor, type ):
		OrderTwoPhenotype.__init__(self, name, k, A, B, max_rate, aff_factor, type)
		self.A = A
		self.B = B
		self.C = C
		self.aff_factor = aff_factor
		self.type = type
		self.max_rate = max_rate
		self.k = k

	def react( self ):
		antigen = self.A.return_antigen()
		agg_rate = self.k * self.A.size()
		base_gene = self.B.select_random_weighted( self.aff_factor, self.max_rate, self.type, agg_rate, antigen )

		self.B.genotype_decrease(base_gene, 1 )
		self.C.genotype_increase(base_gene, 1 )

class MStimulation( OrderTwoPhenotype ):
	def __init__( self, name, k, A, B, C, D, max_rate, aff_factor, type ):
		OrderTwoPhenotype.__init__(self, name, k, A, B, max_rate, aff_factor, type)
		self.A = A
		self.B = B
		self.C = C
		self.D = D    
		self.aff_factor = aff_factor
		self.type = type
		self.max_rate = max_rate
		self.k = k

	def react( self ):
		antigen = self.A.return_antigen()
		agg_rate = self.k * self.A.size()
		base_gene = self.B.select_random_weighted( self.aff_factor, self.max_rate, self.type, agg_rate, antigen )

		self.B.genotype_decrease(base_gene, 1 )
		if (random.random() < 0.975):        
		    self.C.genotype_increase(base_gene, 1 )
		else:
		    self.D.genotype_increase(base_gene, 1 )        
        
#Name: Formation
#Desc: Defines the spontaneous formation of a B cell
#		i.e.	0 -> B
class Formation( OrderZero ):
	def __init__( self, name, k, A ):
		OrderZero.__init__(self, name, k)
		self.A = A

	def react( self ):
		self.A.increase( 1 )

#Name: Decay
#Desc: Define the first order decay of many components, such as antigen, antibodies, plasma cells
#		i.e.	Ag -> 0
class Decay( OrderOne ):
	def __init__( self, name, k, A ):
		self.A = A
		OrderOne.__init__(self, name, k, A)

	def react( self ):
		self.A.decrease( 1 )


#Name: Viral replication
#Desc: Define the first order replication of viral particles
class Replication( OrderOne ):
	def __init__( self, name, k, A ):
		self.A = A
		OrderOne.__init__(self, name, k, A)

	def react( self ):
		self.A.increase( 1 )

#Name: PopulationDecay
#Desc: Defines the first order decay of a population under a carrying capacity, for example GC B cells
#		i.e. B -> 0 
class PopulationDecay:
	def __init__( self, name, r, k, capacity, A, B, C, D):
		self.name = name
		self.r = r			#replication rate
		self.k = k			#minimum decay rate
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.capacity = float(capacity) 	#carrying capacity
		self.__rate = 0

	def rate( self ):
		self.__rate = 0
		population_size = float(self.A.size())
		k_adjust = self.r * (population_size/self.capacity)
		k = max(self.k, k_adjust)

		self.__rate = k * population_size
		return self.__rate

	def react( self ):
		i = random.random()
		if (i <= 0.33):
			if (self.B.size() > 0): 
				self.B.decrease(1)
		elif (i <= 0.66):
			if (self.C.size() > 0):
				self.C.decrease(1)
		else:
			if (self.D.size() > 0):
				self.D.decrease(1)
                
#Name: TPopulationDecay
#Desc: Defines the first order decay of a population under a carrying capacity
class TPopulationDecay:
	def __init__( self, name, r, k, capacity, A, B):
		self.name = name
		self.r = r			#replication rate
		self.k = k			#minimum decay rate
		self.A = A
		self.B = B
		self.capacity = float(capacity) 	#carrying capacity
		self.__rate = 0

	def rate( self ):
		self.__rate = 0
		population_size = float(self.A.size())
		k_adjust = self.r * (population_size/self.capacity)
		k = max(self.k, k_adjust)

		self.__rate = k * population_size
		return self.__rate

	def react( self ):
		if (self.B.size() > 0):
			self.B.decrease(1)

#Name: Differentiation 
#Desc: Defines the first order reaction of differentiation into one of three components, for example 
#      a stimulated B cell can divide into a daughter B cell, memory cell, or plasma cell.
#      The division into a daughter B cell has a probability including a somatic mutation.
#	i.e.	B* -> 2B
#		       -> B + Me
#              -> B + Pl 
class Differentiation( OrderOne ):
	def __init__( self, name, k, A, B, C, D, E, F, max_rate ):
		OrderOne.__init__(self, name, k, A)
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.E = E        
		self.F = F   
		self.max_rate = max_rate

	def react( self ):
		gene_A = self.A.select_random()
		self.A.genotype_decrease(gene_A, 1)

		iso_rate = isotype_rate
		mut_rate = (mutation_rate + mutation_rate/gene_vocab) * (1-lethal_fraction)
		diff_rate = differentiation_rate

		for i in range(2):
			new_gene = gene_A
			r = random.random()
			if ( r <= iso_rate ):
				new_gene = IsotypeSwitch( gene_A )
				self.B.genotype_increase(new_gene, 1)
			elif ( r <= mut_rate+iso_rate ): 
				new_gene = GeneMutate( new_gene )
				self.B.genotype_increase(new_gene, 1)
			elif ( r <= mut_rate+iso_rate+reverse_rate): 
				self.C.genotype_increase(new_gene, 1)
			elif (r <= mut_rate+iso_rate+reverse_rate+diff_rate):
				if(random.random() < 0.2):
					self.D.genotype_increase(gene_A, 1)
				else:
				    s = random.random()
				    if(s < 0.975):
					    self.E.genotype_increase(gene_A, 1)
				    else:
					    self.F.genotype_increase(gene_A, 1)                        
			else:
				self.B.genotype_increase(new_gene, 1)

class BDifferentiation( OrderOne ):
	def __init__( self, name, k, A, B, max_rate ):
		OrderOne.__init__(self, name, k, A)
		self.A = A
		self.B = B      
		self.max_rate = max_rate

	def react( self ):
		gene_A = self.A.select_random()
		self.A.genotype_decrease(gene_A, 1)

		iso_rate = isotype_rate
		mut_rate = (mutation_rate + mutation_rate/gene_vocab) * (1-lethal_fraction)
		diff_rate = differentiation_rate

		for i in range(2):
			new_gene = gene_A
			r = random.random()
			if ( r <= iso_rate ):
				new_gene = IsotypeSwitch( gene_A )
				self.B.genotype_increase(new_gene, 1)
			elif ( r <= mut_rate+iso_rate ): 
				new_gene = GeneMutate( new_gene )
				self.B.genotype_increase(new_gene, 1)                       
			else:
				self.B.genotype_increase(new_gene, 1)
			
#Name: Production
#Desc: Defines the first order production of one component by another component
#      for example plasma cells producing antibodies.
#	i.e. Pl -> Pl + Ab
class Production( OrderOne ):
	def __init__( self, name, k, A, B ):
		OrderOne.__init__(self, name, k, A)
		self.A = A
		self.B = B

	def react( self ):
		base_gene = self.A.select_random()
		self.B.genotype_increase(base_gene, 1)


class LLPCProduction( OrderTwo ):
	def __init__( self, name, k, A, B, C):
		OrderTwo.__init__(self, name, k, A, B)
		self.A = A
		self.B = B
		self.C = C

	def react( self ):
		base_gene = self.B.select_random()
		self.C.genotype_increase(base_gene, 1 )

#Name: Clearance
#Desc: Defines the second order clearance of one component by another component
#      through Ab/BCR binding. For example antibody-based clearance of an antigen
#	i.e. Ab + Ag -> Ab 
class AbClearance( OrderTwoPhenotype ):
	def __init__( self, name, k, A, B, max_rate, aff_factor, type, AbClearcount):
		OrderTwoPhenotype.__init__(self, name, k, A, B, max_rate, aff_factor, type)
		self.A = A
		self.B = B
		self.AbClearcount = AbClearcount
        
	def react( self ):
#		if (self.B.size()/self.A.size() >= 5):
		if (self.B.size() > 900):
			self.A.decrease(1)
#			self.B.decrease(1)
#			base_gene = self.B.select_random()
#			self.B.genotype_decrease(base_gene, 1 )
			self.AbClearcount.increase(1)
#			for i in range(30):           
#				base_gene = self.B.select_random()
#				self.B.genotype_decrease(base_gene, 1 )
        
#Name: Clearance
#Desc: Defines the second order viral clearance by CD8 T-cell
class TClearance( OrderTwo):
	def __init__( self, name, k, A, B, T8Clearcount):
		OrderTwo.__init__(self, name, k, A, B)
		self.A = A
		self.B = B
		self.T8Clearcount = T8Clearcount
        
	def react( self ):
		self.A.decrease(1)
		self.T8Clearcount.increase(1)

        
#Name: Stimulation of B cell by T cell
#Desc: Defines the second order reaction
class TStimulation( OrderTwo ):
	def __init__( self, name, k, A, B, C, D ):
		OrderTwo.__init__(self, name, k, A, B)
		self.A = A
		self.B = B
		self.C = C
		self.D = D
        
	def react( self ):
		base_gene = self.A.select_random()
		self.A.genotype_decrease(base_gene, 1)        
		self.B.decrease(1)    
		self.C.genotype_increase(base_gene, 1)
		self.D.increase(1)

#Name: Stimulation of CD8 T cell by CD4 stimulated T cell
#Desc: Defines the second order reaction
class T8Stimulation( OrderTwo ):
	def __init__( self, name, k, A, B, C):
		OrderTwo.__init__(self, name, k, A, B)
		self.A = A
		self.B = B
		self.C = C
        
	def react( self ):
		self.B.decrease(1)    
		self.C.increase(1)
        
#Name: TCD4 Differentiation by Mutation
#Name: TCD4 Differentiation into memory 
#Desc: Define the first order reaction

class TDifferentiation( OrderOne ):
	def __init__( self, name, k, A, B, C):
		OrderOne.__init__(self, name, k, A)
		self.A = A
		self.B = B
		self.C = C

	def react( self ):

		self.A.decrease(1)        
		r = random.random()
        
		if ( r <= 0.025 ):          
			self.B.increase(2)
		elif ( r <= 0.9 ):          
			self.B.increase(1)       
		else:
			self.C.increase(1)

#Name: TCD4 Differentiation by Mutation
#Name: TCD4 Differentiation into memory 
#Desc: Define the first order reaction

class T8Differentiation( OrderOne ):
	def __init__( self, name, k, A, B):
		OrderOne.__init__(self, name, k, A)
		self.A = A
		self.B = B

	def react( self ):
		self.A.decrease(1)        
		self.B.increase(2)
            
#########################################
## SYSTEM FUNCTIONS 
#########################################
#Name: TotalReaction
#Desc: Contains the entire set of immune reactions that define the system. Calculates
#      the reaction rate and Monte Carlo time step based on the Gillespie algorithm
class TotalReaction:
    def __init__( self ):
        self.reaction_list = []
        self.n = 0
        self.__rate = []
        self.__total_rate = float(0.0)
        
    def __len__(self):
        return len(self.reaction_list)       

    def add_reaction( self, reaction ):
        self.reaction_list.append( reaction )
        self.n = len(self.reaction_list)

    def remove_reaction( self, reaction ):
        self.reaction_list.remove( reaction )
        self.n = len(self.reaction_list)
        
    def rate( self ):
        self.__rate = []
        for i in range(0,self.n):
            self.__rate.append( self.reaction_list[i].rate() )
        self.__total_rate = sum( self.__rate )
        
    def MC_TimeStep( self ):
        r = random.random()
        self.rate()
        if (self.__total_rate > 0.0):
            dt = (1/self.__total_rate)*math.log(1/r)
            return dt        
        else:
            return 0.0
            
    def MC_React( self ):
        r = random.random()
#       self.rate()

        outcome = float( 0.0 )
        
        for i in range(0,self.n):
            if (self.__total_rate > 0.0): 
                outcome += (self.__rate[i]/self.__total_rate)
            else: 
                outcome += 0.0
            if (r <= outcome):
                self.reaction_list[i].react()
                break

#Name: AntigenFileInput
#Desc: Inputs the antigen information in terms of antigen and epitope parameters
#      Returns an array of antigens
def ReadAgFile( filename ):
    f = open( filename, 'r')
    f.readline()

    ag_index = -1
    antigen_list = []
    old_ag = "na"
    
    for line in f:
        d1 = line.rsplit('\n')  #strip end of line
        d2 = d1[0].split(';')   #semi-colon delimited
        ag_name = d2[0]
        
        if (old_ag != ag_name):
            old_ag = ag_name
            ag_index = ag_index + 1
            antigen_list.append( AntigenType(ag_name, ag_index))
			
        ep_name = str(d2[1])
        ep_seq = str(d2[2])
        ep_immunogenicity = float(d2[3])
        ep_clearance = float(d2[4])
        antigen_list[ ag_index ].add_epitope( Epitope(ep_name, ep_seq, ep_immunogenicity, ep_clearance))
        
    return antigen_list

#Name: FileOutput 
#Desc: Outputs the simulation data into a data file, at defined time intervals. The simulation 
#      data is mainly the individual populations of each component of the system, as well
#      as subpopulations broken down by phenotype. 
class FileOutput:
	def __init__( self, data_file, step, population, antigen_list ):
		self.time = float( 0.0 )
		self.data_file = data_file
		self.step = step
		self.population = population
		self._antigen_list = antigen_list
		self._antigen1 = self._antigen_list[0]
		self._antigen2 = self._antigen_list[1]
		self._antigen3 = self._antigen_list[2]
		self._antigen4 = self._antigen_list[3]        

	def write( self, time ):
		if (time == 0.0 or time > (self.time + self.step)):
			line_out = ''
            #  tau = 8.0, check!
			line_out += str( time*8.0/24.0 ) +';'

			#for i in range(0,len(self.population)):
			#	line_out += str(self.population[i].size())+';'

#			line_out += 'Clearance;'
			line_out += str(self.population[17].size())+';' #AbClear
			line_out += str(self.population[18].size())+';' #T8Clear
			
#			line_out += 'V;'
			for i in range(6, 10):
				line_out += str(self.population[i].size())+';'
#			line_out += 'T;'
			line_out += str(self.population[10].size()+self.population[11].size()+self.population[12].size())+';' #Total Tcd4            
			line_out += str(self.population[10].size())+';' #Tcd4
			line_out += str(self.population[11].size())+';' #Tcd4_stim
			line_out += str(self.population[12].size())+';' #Tcd4_me
			line_out += str(self.population[13].size()+self.population[16].size())+';' #Total Tcd8          
			line_out += str(self.population[13].size())+';' #Tcd8
			line_out += str(self.population[16].size())+';' #Tcd8_stim
            
#			line_out += 'B;'
			line_out += str(self.population[0].size()+self.population[1].size()+self.population[14].size())+';' #Total B            
			line_out += str(self.population[0].size())+';' #GC
			line_out += str(self.population[1].size())+';' #Antigen_stim_B
			line_out += str(self.population[2].size())+';' #Memory
			line_out += str(self.population[3].size())+';' #Plasma SL
			line_out += str(self.population[15].size())+';' #Plasma LL            
			line_out += str(self.population[4].size())+';' #Antibody
			line_out += str(self.population[5].size())+';' #Naive
			line_out += str(self.population[14].size())+';'#Tcell_stim_B			

#			line_out += 'Bphen;'
			for i in range(0,8):
				line_out += str(self.population[0].phenotype_size(i, self._antigen1.ID()))+';'	
			for i in range(0,8):
				line_out += str(self.population[1].phenotype_size(i, self._antigen1.ID()))+';'
			for i in range(0,8):
				line_out += str(self.population[14].phenotype_size(i, self._antigen1.ID()))+';'                
			for i in range(0,8):
				line_out += str(self.population[2].phenotype_size(i, self._antigen1.ID()))+';'
			for i in range(0,8):
				output = []
				output.append(self.population[3].phenotype_size(i, self._antigen1.ID()))
				output.append(self.population[15].phenotype_size(i, self._antigen1.ID()))
				line_out += str(max(output))+';'
                
#			line_out += 'ABphen;'
			for i in range(0,8):
				line_out += str(self.population[4].phenotype_size(i, self._antigen1.ID()))+';'	
                
#			line_out += 'ABcross_specificity;'
			output = self.population[4].calc_crossreactivity_specificity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'

#			line_out += 'ABcross;'
			output = self.population[4].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'            
                      
#			line_out += 'Ncross;'
			output = self.population[5].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
            
#			line_out += 'GCcross;'
			output = self.population[0].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
            
#			line_out += 'B*cross;'
			output = self.population[1].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
            
#			line_out += 'B**cross;'
			output = self.population[14].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
			
#			line_out += 'Mcross;'
			output = self.population[2].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
            
#			line_out += 'Pcross;'
			output = self.population[15].calc_crossreactivity(self._antigen_list)
			for m in range(0,len(output)):
				line_out += str(output[m])+';'
            
#			line_out += 'ABtrans;'
			output = self.population[4].calc_transcend(self._antigen_list)
			#for i in range(0, len(self._antigen_list)+1):
			for i in range(0, len(self._antigen_list)):
				line_out += str(output[i])+';'

#			line_out += 'Btrans;'		
			output0 = self.population[0].calc_transcend(self._antigen_list)
			output1 = self.population[1].calc_transcend(self._antigen_list)
			output2 = self.population[2].calc_transcend(self._antigen_list)
			output3 = self.population[14].calc_transcend(self._antigen_list)
			#for i in range(0, len(self._antigen_list)+1):
			for i in range(0, len(self._antigen_list)):
				output = [] 
				output.append(output0[i])
				output.append(output1[i])
				output.append(output2[i])
				output.append(output3[i])
				line_out += str(max(output))+';'

#			line_out += 'Bdiv;'
			output = []
			n = 0.25
			output.append(self.population[0].diversity2(n))
			output.append(self.population[1].diversity2(n))
			output.append(self.population[2].diversity2(n))
			output.append(self.population[14].diversity2(n))
			line_out += str(max(output))+';'

			n = 0.50
			output.append(self.population[0].diversity2(n))
			output.append(self.population[1].diversity2(n))
			output.append(self.population[2].diversity2(n))
			output.append(self.population[14].diversity2(n))
			line_out += str(max(output))+';'

			n = 0.75
			output.append(self.population[0].diversity2(n))
			output.append(self.population[1].diversity2(n))
			output.append(self.population[2].diversity2(n))
			output.append(self.population[14].diversity2(n))
			line_out += str(max(output))+';'
			
			n = 1.0
			output.append(self.population[0].diversity2(n))
			output.append(self.population[1].diversity2(n))
			output.append(self.population[2].diversity2(n))
			output.append(self.population[14].diversity2(n))
			line_out += str(max(output))+';'

#			line_out += 'AB_ep1A_phen;'
			output = self.population[4].calc_cross(self._antigen1, 0)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep1B_phen;'
			output = self.population[4].calc_cross(self._antigen1, 1)	
			for i in range(0,4):
				line_out += str(output[i])+';'
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep1C_phen;'
			output = self.population[4].calc_cross(self._antigen1, 2)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep1D_phen;'
			output = self.population[4].calc_cross(self._antigen1, 3)	
			for i in range(0,4):
				line_out += str(output[i])+';'
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep2A_phen;'
			output = self.population[4].calc_cross(self._antigen2, 0)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
			
#			line_out += 'AB_ep2B_phen;'
			output = self.population[4].calc_cross(self._antigen2, 1)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep2C_phen;'
			output = self.population[4].calc_cross(self._antigen2, 2)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
			
#			line_out += 'AB_ep2D_phen;'
			output = self.population[4].calc_cross(self._antigen2, 3)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep3A_phen;'
			output = self.population[4].calc_cross(self._antigen3, 0)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep3B_phen;'
			output = self.population[4].calc_cross(self._antigen3, 1)	
			for i in range(0,4):
				line_out += str(output[i])+';'
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep3C_phen;'
			output = self.population[4].calc_cross(self._antigen3, 2)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep3D_phen;'
			output = self.population[4].calc_cross(self._antigen3, 3)	
			for i in range(0,4):
				line_out += str(output[i])+';'
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep4A_phen;'
			output = self.population[4].calc_cross(self._antigen4, 0)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
			
#			line_out += 'AB_ep4B_phen;'
			output = self.population[4].calc_cross(self._antigen4, 1)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
            
#			line_out += 'AB_ep4C_phen;'
			output = self.population[4].calc_cross(self._antigen4, 2)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'
			
#			line_out += 'AB_ep4D_phen;'
			output = self.population[4].calc_cross(self._antigen4, 3)	
			for i in range(0,4):
				line_out += str(output[i])+';'	
			line_out += str(output[0]+output[1]+output[2]+output[3])+';'            
            
#			line_out += 'AB_Ig_MorG;'
			output = self.population[4].calc_isotype()	
			for i in range(0,2):
				line_out += str(output[i])+';'	
			
#			print(line_out)
			print('Time ', time/3.0)
			line_out += '\n'
			self.f.write( line_out )
			self.time = time

	def start( self ):
		self.f = open( self.data_file, 'w' )
		line_out = ''
		line_out += "Time;"

#		line_out += 'Clearance;'
		line_out += str(self.population[17].name())+';'
		line_out += str(self.population[18].name())+';'
        
#		line_out += 'V;'
		for i in range(6, 10):
				line_out += str(self.population[i].name())+';'
                
#		line_out += 'T;'
		line_out += 'Total Tcd4;'
		line_out += str(self.population[10].name())+';'
		line_out += str(self.population[11].name())+';'
		line_out += str(self.population[12].name())+';'
		line_out += 'Total Tcd8;'         
		line_out += str(self.population[13].name())+';'        
		line_out += str(self.population[16].name())+';'
        
#		line_out += 'B;'
		line_out += 'Total B;'        
		line_out += str(self.population[0].name())+';'
		line_out += str(self.population[1].name())+';'
		line_out += str(self.population[2].name())+';'
		line_out += str(self.population[3].name())+';'
		line_out += str(self.population[15].name())+';'       
		line_out += str(self.population[4].name())+';'
		line_out += str(self.population[5].name())+';'
		line_out += str(self.population[14].name())+';'
        
#		line_out += 'Bphen;'
		tag = "B_"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'
		tag = "Bstim_"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'
		tag = "Bact_"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'
		tag = "Bmem_"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'
		tag = "Bplas_"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'       
#		line_out += 'ABphen;'
		tag = "AB"
		for i in range(0,8):
			line_out += tag+"_"+str(i)+';'
            
#		line_out += 'ABcross_specificity;'
		tag = "AB_cross1234_1"
		line_out += tag+';'
		tag = "AB_cross1234_2"
		line_out += tag+';'
		tag = "AB_cross1234_3"
		line_out += tag+';'
		tag = "AB_cross1234_4"
		line_out += tag+';'        
		tag = "AB_cross123_1"
		line_out += tag+';'        
		tag = "AB_cross123_2"
		line_out += tag+';'        
		tag = "AB_cross123_3"
		line_out += tag+';'        
		tag = "AB_cross123_4"        
		line_out += tag+';'
		tag = "AB_cross124_1"
		line_out += tag+';'
		tag = "AB_cross124_2"
		line_out += tag+';'
		tag = "AB_cross124_3"
		line_out += tag+';'
		tag = "AB_cross124_4"
		line_out += tag+';'       
		tag = "AB_cross134_1"
		line_out += tag+';'
		tag = "AB_cross134_2"
		line_out += tag+';' 
		tag = "AB_cross134_3"
		line_out += tag+';' 
		tag = "AB_cross134_4"
		line_out += tag+';'        
		tag = "AB_cross234_1"
		line_out += tag+';'
		tag = "AB_cross234_2"
		line_out += tag+';'
		tag = "AB_cross234_3"
		line_out += tag+';'
		tag = "AB_cross234_4"
		line_out += tag+';'        
		tag = "AB_cross12_1"
		line_out += tag+';'
		tag = "AB_cross12_2"
		line_out += tag+';'
		tag = "AB_cross12_3"
		line_out += tag+';'
		tag = "AB_cross12_4"
		line_out += tag+';'        
		tag = "AB_cross13_1"
		line_out += tag+';'		
		tag = "AB_cross13_2"
		line_out += tag+';'
		tag = "AB_cross13_3"
		line_out += tag+';' 
		tag = "AB_cross13_4"
		line_out += tag+';'        
		tag = "AB_cross14_1"
		line_out += tag+';'
		tag = "AB_cross14_2"
		line_out += tag+';'
		tag = "AB_cross14_3"
		line_out += tag+';'
		tag = "AB_cross14_4"
		line_out += tag+';'        
		tag = "AB_cross23_1"
		line_out += tag+';'
		tag = "AB_cross23_2"
		line_out += tag+';' 
		tag = "AB_cross23_3"
		line_out += tag+';' 
		tag = "AB_cross23_4"
		line_out += tag+';'         
		tag = "AB_cross24_1"
		line_out += tag+';'
		tag = "AB_cross24_2"
		line_out += tag+';'
		tag = "AB_cross24_3"
		line_out += tag+';'
		tag = "AB_cross24_4"
		line_out += tag+';'        
		tag = "AB_cross34_1"
		line_out += tag+';'
		tag = "AB_cross34_2"
		line_out += tag+';'
		tag = "AB_cross34_3"
		line_out += tag+';'
		tag = "AB_cross34_4"
		line_out += tag+';'        
		tag = "AB_spec1_1"
		line_out += tag+';'
		tag = "AB_spec1_2"
		line_out += tag+';'
		tag = "AB_spec1_3"
		line_out += tag+';'
		tag = "AB_spec1_4"
		line_out += tag+';'        
		tag = "AB_spec2_1"
		line_out += tag+';'
		tag = "AB_spec2_2"
		line_out += tag+';'
		tag = "AB_spec2_3"
		line_out += tag+';'
		tag = "AB_spec2_4"
		line_out += tag+';'        
		tag = "AB_spec3_1"
		line_out += tag+';'
		tag = "AB_spec3_2"
		line_out += tag+';'
		tag = "AB_spec3_3"
		line_out += tag+';'
		tag = "AB_spec3_4"
		line_out += tag+';'        
		tag = "AB_spec4_1"
		line_out += tag+';'
		tag = "AB_spec4_2"
		line_out += tag+';'
		tag = "AB_spec4_3"
		line_out += tag+';'
		tag = "AB_spec4_4"
		line_out += tag+';'
        
#		line_out += 'ABcross;'
		tag = "AB_cross1234"
		line_out += tag+';'
		tag = "AB_cross123"
		line_out += tag+';'
		tag = "AB_cross124"
		line_out += tag+';'
		tag = "AB_cross134"
		line_out += tag+';'        
		tag = "AB_cross234"
		line_out += tag+';'
		tag = "AB_cross12"
		line_out += tag+';'
		tag = "AB_cross13"
		line_out += tag+';'
		tag = "AB_cross14"
		line_out += tag+';'
		tag = "AB_cross23"
		line_out += tag+';' 
		tag = "AB_cross24"
		line_out += tag+';'
		tag = "AB_cross34"
		line_out += tag+';'
		tag = "AB_spec1"
		line_out += tag+';'
		tag = "AB_spec2"
		line_out += tag+';'
		tag = "AB_spec3"
		line_out += tag+';'
		tag = "AB_spec4"
		line_out += tag+';'
		tag = "AB_other"
		line_out += tag+';'
        
#		line_out += 'Ncross;'
		tag = "N_cross1234"
		line_out += tag+';'
		tag = "N_cross123"
		line_out += tag+';'
		tag = "N_cross124"
		line_out += tag+';'
		tag = "N_cross134"
		line_out += tag+';'        
		tag = "N_cross234"
		line_out += tag+';'
		tag = "N_cross12"
		line_out += tag+';'
		tag = "N_cross13"
		line_out += tag+';'
		tag = "N_cross14"
		line_out += tag+';'
		tag = "N_cross23"
		line_out += tag+';' 
		tag = "N_cross24"
		line_out += tag+';'
		tag = "N_cross34"
		line_out += tag+';'
		tag = "N_spec1"
		line_out += tag+';'
		tag = "N_spec2"
		line_out += tag+';'
		tag = "N_spec3"
		line_out += tag+';'
		tag = "N_spec4"
		line_out += tag+';'
		tag = "N_other"
		line_out += tag+';'

#		line_out += 'GCcross;'
		tag = "GC_cross1234"
		line_out += tag+';'
		tag = "GC_cross123"
		line_out += tag+';'
		tag = "GC_cross124"
		line_out += tag+';'
		tag = "GC_cross134"
		line_out += tag+';'        
		tag = "GC_cross234"
		line_out += tag+';'
		tag = "GC_cross12"
		line_out += tag+';'
		tag = "GC_cross13"
		line_out += tag+';'
		tag = "GC_cross14"
		line_out += tag+';'
		tag = "GC_cross23"
		line_out += tag+';' 
		tag = "GC_cross24"
		line_out += tag+';'
		tag = "GC_cross34"
		line_out += tag+';'
		tag = "GC_spec1"
		line_out += tag+';'
		tag = "GC_spec2"
		line_out += tag+';'
		tag = "GC_spec3"
		line_out += tag+';'
		tag = "GC_spec4"
		line_out += tag+';'
		tag = "GC_other"
		line_out += tag+';'
 
#		line_out += 'B*cross;'
		tag = "B*_cross1234"
		line_out += tag+';'
		tag = "B*_cross123"
		line_out += tag+';'
		tag = "B*_cross124"
		line_out += tag+';'
		tag = "B*_cross134"
		line_out += tag+';'        
		tag = "B*_cross234"
		line_out += tag+';'
		tag = "B*_cross12"
		line_out += tag+';'
		tag = "B*_cross13"
		line_out += tag+';'
		tag = "B*_cross14"
		line_out += tag+';'
		tag = "B*_cross23"
		line_out += tag+';' 
		tag = "B*_cross24"
		line_out += tag+';'
		tag = "B*_cross34"
		line_out += tag+';'
		tag = "B*_spec1"
		line_out += tag+';'
		tag = "B*_spec2"
		line_out += tag+';'
		tag = "B*_spec3"
		line_out += tag+';'
		tag = "B*_spec4"
		line_out += tag+';'
		tag = "B*_other"
		line_out += tag+';'

#		line_out += 'B**cross;'
		tag = "B**_cross1234"
		line_out += tag+';'
		tag = "B**_cross123"
		line_out += tag+';'
		tag = "B**_cross124"
		line_out += tag+';'
		tag = "B**_cross134"
		line_out += tag+';'        
		tag = "B**_cross234"
		line_out += tag+';'
		tag = "B**_cross12"
		line_out += tag+';'
		tag = "B**_cross13"
		line_out += tag+';'
		tag = "B**_cross14"
		line_out += tag+';'
		tag = "B**_cross23"
		line_out += tag+';' 
		tag = "B**_cross24"
		line_out += tag+';'
		tag = "B**_cross34"
		line_out += tag+';'
		tag = "B**_spec1"
		line_out += tag+';'
		tag = "B**_spec2"
		line_out += tag+';'
		tag = "B**_spec3"
		line_out += tag+';'
		tag = "B**_spec4"
		line_out += tag+';'
		tag = "B**_other"
		line_out += tag+';'
        
#		line_out += 'Mcross;'
		tag = "M_cross1234"
		line_out += tag+';'
		tag = "M_cross123"
		line_out += tag+';'
		tag = "M_cross124"
		line_out += tag+';'
		tag = "M_cross134"
		line_out += tag+';'        
		tag = "M_cross234"
		line_out += tag+';'
		tag = "M_cross12"
		line_out += tag+';'
		tag = "M_cross13"
		line_out += tag+';'
		tag = "M_cross14"
		line_out += tag+';'
		tag = "M_cross23"
		line_out += tag+';' 
		tag = "M_cross24"
		line_out += tag+';'
		tag = "M_cross34"
		line_out += tag+';'
		tag = "M_spec1"
		line_out += tag+';'
		tag = "M_spec2"
		line_out += tag+';'
		tag = "M_spec3"
		line_out += tag+';'
		tag = "M_spec4"
		line_out += tag+';'
		tag = "M_other"
		line_out += tag+';'
				
#		line_out += 'Pcross;'
		tag = "P_cross1234"
		line_out += tag+';'
		tag = "P_cross123"
		line_out += tag+';'
		tag = "P_cross124"
		line_out += tag+';'
		tag = "P_cross134"
		line_out += tag+';'        
		tag = "P_cross234"
		line_out += tag+';'
		tag = "P_cross12"
		line_out += tag+';'
		tag = "P_cross13"
		line_out += tag+';'
		tag = "P_cross14"
		line_out += tag+';'
		tag = "P_cross23"
		line_out += tag+';' 
		tag = "P_cross24"
		line_out += tag+';'
		tag = "P_cross34"
		line_out += tag+';'
		tag = "P_spec1"
		line_out += tag+';'
		tag = "P_spec2"
		line_out += tag+';'
		tag = "P_spec3"
		line_out += tag+';'
		tag = "P_spec4"
		line_out += tag+';'
		tag = "P_other"
		line_out += tag+';'
        
#		line_out += 'ABtrans;'
		tag = "AB_trans"
		for i in range(0, len(self._antigen_list)):
			line_out += tag+"_"+str(i)+';'

#		line_out += 'Btrans;'
		tag = "B_trans"
		for i in range(0, len(self._antigen_list)):
			line_out += tag+"_"+str(i)+';'
		
#		line_out += 'Bdiv;'
		tag = "Bdiv_25"
		line_out += tag+';'
		tag = "Bdiv_50"
		line_out += tag+';'
		tag = "Bdiv_75"
		line_out += tag+';'
		tag = "Bdiv_100"
		line_out += tag+';'
		
#		line_out += 'AB_ep1A_phen;'
		tag = "ABE1A"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep1Aphen;'
        
#		line_out += 'AB_ep1B_phen;'
		tag = "ABE1B"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep1Bphen;'

#		line_out += 'AB_ep1C_phen;'
		tag = "ABE1C"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep1Cphen;'

#		line_out += 'AB_ep1D_phen;'
		tag = "ABE1D"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep1Dphen;'
        
#		line_out += 'AB_ep2A_phen;'
		tag = "ABE2A"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep2Aphen;'
        
#		line_out += 'AB_ep2B_phen;'
		tag = "ABE2B"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep2Bphen;'
 
#		line_out += 'AB_ep2C_phen;'
		tag = "ABE2C"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep2Cphen;'
        
#		line_out += 'AB_ep2D_phen;'
		tag = "ABE2D"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep2Dphen;'
 
#		line_out += 'AB_ep3A_phen;'
		tag = "ABE3A"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep3Aphen;'
        
#		line_out += 'AB_ep3B_phen;'
		tag = "ABE3B"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep3Bphen;'

#		line_out += 'AB_ep3C_phen;'
		tag = "ABE3C"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep3Cphen;'

#		line_out += 'AB_ep3D_phen;'
		tag = "ABE3D"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep3Dphen;'
        
#		line_out += 'AB_ep4A_phen;'
		tag = "ABE4A"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep4Aphen;'
        
#		line_out += 'AB_ep4B_phen;'
		tag = "ABE4B"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep4Bphen;'
 
#		line_out += 'AB_ep4C_phen;'
		tag = "ABE4C"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep4Cphen;'
        
#		line_out += 'AB_ep4D_phen;'
		tag = "ABE4D"
		for i in range(0,4):
			line_out += tag+"_"+str(i)+';'
		line_out += 'ABep4Dphen;'
        
#		line_out += 'AB_Ig_MorG;'
		line_out += 'AB_IgM;'
		line_out += 'AB_IgG;'           
		
		line_out += '\n'
		self.f.write( line_out )	#outputs header
		self.write( 0.0 )		#writes starting conditions

	def finish( self ):
		self.f.close()
		
	def set_time( t ):
		self.time = t


	def gene_out( self ):
		gene_file = self.data_file + ".gen"
		self.f = open( gene_file, 'w' )
	
		gene_list = self.population[4].return_gene_list() 
		size_list = self.population[4].return_size_list() 
		for i in range(0, len(gene_list)):
			line_out = ''
			line_out += str(gene_list[i]) + ';'
			line_out += str(GeneEpitope(gene_list[i], self._antigen_list))+';'
			line_out += str(size_list[i])+';'
			line_out += '\n'
#			print(line_out)
			self.f.write( line_out )	#outputs header
		self.f.close()
