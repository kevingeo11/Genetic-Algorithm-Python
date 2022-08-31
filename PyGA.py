import numpy as np
from numpy.random import randint
from random import random as rnd

def get_gene(ll, ul, resolution):
	if resolution == 0:
		gene = int(round(rnd()*(ul-ll)+ll, resolution))
	else:
		gene = round(rnd()*(ul-ll)+ll, resolution)
	return gene

def get_chromosome(n_genes, ll, ul, resolution):
	chromosome = [get_gene(ll[x], ul[x], resolution[x]) for x in range(n_genes)]
	return chromosome 

def get_population(n_chromosomes, n_genes, ll, ul, resolution):
	return [get_chromosome(n_genes, ll, ul, resolution) for x in range(n_chromosomes)]

def create_first_generation(pop, fitness_fun):
	fitness = [fitness_fun(pop[x]) for x in range(len(pop))]
	sorted_population_fitness = sorted([[pop[x], fitness[x]]for x in range(len(pop))], key=lambda x: x[1])
	population = [sorted_population_fitness[x][0] for x in range(len(sorted_population_fitness))]
	fitness = [sorted_population_fitness[x][1] for x in range(len(sorted_population_fitness))]
	return {'Chromosomes': population, 'Fitness': fitness}

def fitness_saturation(max_fitness, number_of_similarity):
	result = False
	similarity = 0
	if len(max_fitness) > number_of_similarity:
		for i in range(1,number_of_similarity):
			if max_fitness[len(max_fitness) - i] == max_fitness[len(max_fitness) - i - 1]:
				similarity += 1
			else:
				similarity = 0
		if similarity == number_of_similarity-1:
			result = True
		return result
	else:
		return result

def selection(generation, method='Fittest Half'):
	
	generation['Normalized Fitness'] = sorted([generation['Fitness'][x]/sum(generation['Fitness']) for x in range(len(generation['Fitness']))], reverse = True)
	'''
	Normalized Fitness in decending order
	But chromosomes are in ascending order according to fitness
	'''
	generation['Cumulative Sum'] = np.array(generation['Normalized Fitness']).cumsum()

	if method == 'Roulette Wheel':
		selected = []
		for x in range(len(generation['Individuals'])//2):
			selected.append(roulette(generation
				['Cumulative Sum'], rnd()))
			while len(set(selected)) != len(selected):
				selected[x] = \
					(roulette(generation['Cumulative Sum'], rnd()))
		selected = {'Individuals': 
			[generation['Individuals'][int(selected[x])]
				for x in range(len(generation['Individuals'])//2)]
				,'Fitness': [generation['Fitness'][int(selected[x])]
				for x in range(
					len(generation['Individuals'])//2)]}

	elif method == 'Fittest Half':
		'''
		Takes the last half of the chromosomes i.e. half having the most fitness
		'''
		selected_chromosomes = [generation['Chromosomes'][-x-1] for x in range(len(generation['Chromosomes'])//2)]
		selected_fitness = [generation['Fitness'][-x-1] for x in range(len(generation['Chromosomes'])//2)]
		selected = {'Chromosomes': selected_chromosomes, 'Fitness': selected_fitness}

	return selected

def crossover(selected, method = 'Fittest'):
	'''
	Creating parents in mating pool for generating offspring
	Chromosomes come in decending order of fitness
	'''
	Chromosomes = selected['Chromosomes']
	fitness = selected['Fitness']
	if method == 'Fittest':
		parents = [[Chromosomes[x], Chromosomes[x+1]] for x in range(len(Chromosomes)//2 if len(Chromosomes)%2==0 else len(Chromosomes)//2 + 1)]
	if method == 'Random':
		parents = []
		for x in range(len(Chromosomes)//2):
			parents.append(
				[Chromosomes[randint(0,(len(Chromosomes)-1))],
				 Chromosomes[randint(0,(len(Chromosomes)-1))]])
			while parents[x][0] == parents[x][1]:
				parents[x][1] = Chromosomes[
					randint(0,(len(Chromosomes)-1))]
	'''
	Cross over begins
	'''
	childrens = []
	for parent in parents:
		if len(parent[0]) == 1:
			offspring = parent
		else:
			pivot_point = randint(1, len(parent[0]))
			offspring = [parent[0][0:pivot_point] + parent[1][pivot_point:]]
			offspring.append(parent[1][0:pivot_point] + parent[0][pivot_point:])
		
		childrens.append(offspring[0])
		childrens.append(offspring[1])
	offsprings = selected['Chromosomes'] + childrens

	return offsprings

def mutation(offsprings, length, ll, ul, resolution):

	mutated_chromosomes = []
	for i in range(length):
		mutated_chromosome = offsprings[i].copy()
		if randint(0,2) == 0:
			if len(mutated_chromosome) == 1:
				mutated_chromosome = offsprings[i].copy()
			else:
				if len(mutated_chromosome) > 3:
					muatation_rate = min(3, len(mutated_chromosome)//2)
				else:
					muatation_rate = 1 
				
				for x in range(muatation_rate):
					r = randint(0, len(mutated_chromosome))
					mutated_chromosome[r] = get_gene(ll[r], ul[r], resolution[r])
		
		mutated_chromosomes.append(mutated_chromosome)

	return mutated_chromosomes

def next_generation(gen, fitness_fun, ll, ul, resolution, elit_flag):
	if elit_flag == True:
		elit = {}
		elit_selected = {}
		next_generation = {}
		elit['Chromosomes'] = gen['Chromosomes'].pop(-1)
		elit['Fitness'] = gen['Fitness'].pop(-1)
		selected = selection(gen)
		elit_selected['Chromosomes'] = [elit['Chromosomes']] + selected['Chromosomes']
		elit_selected['Fitness'] = [elit['Fitness']] + selected['Fitness']
		offsprings = crossover(elit_selected)

		if len(offsprings) < len(gen['Chromosomes']):
			print('Error 006 : Size mismatch in offsprings')
			exit()

		mutated_offsprings = mutation(offsprings, len(gen['Chromosomes']), ll, ul, resolution)
		new_generation = mutated_offsprings + [elit['Chromosomes']]
		new_generation_fitness = [fitness_fun(new_generation[x]) for x in range(len(new_generation))]
		sorted_next_gen = sorted([[new_generation[x], new_generation_fitness[x]] for x in range(len(new_generation))], key=lambda x: x[1])
		next_generation['Chromosomes'] = [sorted_next_gen[x][0] for x in range(len(sorted_next_gen))]
		next_generation['Fitness'] = [sorted_next_gen[x][1] for x in range(len(sorted_next_gen))]
		'''
		appending back elite chromosome
		'''
		gen['Chromosomes'].append(elit['Chromosomes'])
		gen['Fitness'].append(elit['Fitness'])
	else:
		next_generation = {}
		selected = selection(gen)
		offsprings = crossover(selected)

		if len(offsprings) < len(gen['Chromosomes']):
			print('Error 005 : Size mismatch in offsprings')
			exit()

		mutated_offsprings = mutation(offsprings, len(gen['Chromosomes']), ll, ul, resolution)
		new_generation = mutated_offsprings
		new_generation_fitness = [fitness_fun(new_generation[x]) for x in range(len(new_generation))]
		sorted_next_gen = sorted([[new_generation[x], new_generation_fitness[x]] for x in range(len(new_generation))], key=lambda x: x[1])
		next_generation['Chromosomes'] = [sorted_next_gen[x][0] for x in range(len(sorted_next_gen))]
		next_generation['Fitness'] = [sorted_next_gen[x][1] for x in range(len(sorted_next_gen))]
		
	return next_generation

def GeneticAlgorithm(fitness_fun, n_par, ll, ul, n_pop = 200, n_generations = 100, n_saturated_generations = 50, resolution = None, elit_flag = False):
	'''
	Creating the 1st Population
	'''
	print(' ')
	if n_par != len(ll):
		print("Error 001 : Invalid length no of parameter !~ lower limit length")
		exit()
	if n_par != len(ul):
		print("Error 002 : Invalid length no of parameter !~ upper limit length")
		exit()
	if n_pop < 3:
		print("Error 003 : Invalid size of population")
		exit()
	if resolution is None:
		resolution = [2 for x in range(n_par)]
	else:
		if len(resolution) != n_par:
			print("Error 004 : Invalid resolution length")
			exit()
		for i in range(n_par):
			if resolution[i] < 0:
				print("Error 005 : Invalid resolution type")
				exit()


	gen = []
	fitness_max = []
	pop = get_population(n_pop, n_par, ll, ul, resolution)
	gen.append(create_first_generation(pop, fitness_fun))
	fitness_max.append(max(gen[0]['Fitness']))
	gen_count = 1
	flag = True
	for gen_count in range(1, n_generations):
		print('Gen {} : best fitness = {}'.format(gen_count,max(fitness_max)))
		print('Best Chromosome : ', gen[-1]['Chromosomes'][-1])
		if fitness_saturation(fitness_max, n_saturated_generations) == True:
			print("*** Saturation Acheived ***")
			print('Population = ', gen[-1]['Chromosomes'])
			print('Best Chromosome : ', gen[-1]['Chromosomes'][-1])
			flag = False
			break
		gen.append(next_generation(gen[-1], fitness_fun, ll, ul, resolution, elit_flag))
		fitness_max.append(max(gen[-1]['Fitness']))
	if flag == True:
		print('Population = ', gen[-1]['Chromosomes'])
		print('Best Chromosome : ', gen[-1]['Chromosomes'][-1])
