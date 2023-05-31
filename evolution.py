from deap import base, creator, tools
import traits_cleaning as traits
import functions 
import seaborn as sns
import matplotlib.pyplot as plt
import random
import multiprocessing
import numpy as np

POLYGENICITY = 'high'
REPLICATES = 50
POP = 1000
GEN = 500
MUTPB =  0.0000168 # This is definitely high
CXPB = 0.5

# Start from a nice straightforward DEAP default

# Create individual type
creator.create("Fitness", base.Fitness, weights = (1.0,)) # PRS is analogous to fitness here
creator.create("PRS", float)
creator.create("Individual", list, fitness=creator.Fitness, PRS = creator.PRS)

toolbox = base.Toolbox()

# Register operators 
toolbox.register("mate", functions.cxDiploid, indpb = 0.5)
toolbox.register("mutate",functions.rareMutation)
toolbox.register("select", functions.selRoulette)
# toolbox.register("select", tools.selRoulette)
# toolbox.register("evaluate",functions.PRS_score, trait_table)


if __name__ == "__main__":
    cpu_count = multiprocessing.cpu_count()
    print("Running on {0} cpu".format(cpu_count))
    pool = multiprocessing.Pool(cpu_count)
    toolbox.register("map",pool.map)
    
    
    retention = []
    starting_fitness_distributions = []
    
    fig, ax = plt.subplots(2, sharex='col')
    ax[0].set(title="PRS")
    ax[1].set(title="Heterozygosity")
    
    for x in range(REPLICATES):
        mean_heights = []
        mean_fits = []
        heterozygosity = []
        if POLYGENICITY == 'low':
            trait_table = traits.traits.sample(50).reset_index(drop=True)
        if POLYGENICITY == 'med':
            trait_table = traits.traits.sample(300).reset_index(drop=True)
        if POLYGENICITY == 'high':
            trait_table = traits.traits.sample(2000).reset_index(drop=True)
        if POLYGENICITY == 'all':
            trait_table = traits.traits.sample(2639).reset_index(drop=True)
        
        risk_allele_freqs = list(trait_table['RAF'])
        effect_sizes = list(trait_table['Effect'])
            
        
        toolbox.register("individual", functions.make_individual, creator.Individual, risk_allele_freqs)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate",functions.calc_PGS, effect_sizes)
        
        print("Iteration: ", x)
        pop = toolbox.population(n=POP)
        
        # Evaluate population 
        pgs = list(toolbox.map(toolbox.evaluate, pop))
        # Get mean and standard dev
        starting_pop_mean = np.mean(pgs)
        starting_pop_dev = np.std(pgs)
        
        # Convert to disease risk 
        fitnesses = functions.calc_fitness(pgs, starting_pop_mean, starting_pop_dev)
        starting_fitness_distributions.append(fitnesses)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit,
        
        
        g = 0    
        mutation_retained = True
        
        # fitnesses.sort()
        # plt.hist(fitnesses)
        # plt.title("Starting PRS Distribution") 
        # plt.show()
        # plt.close()
        while g <= GEN and mutation_retained == True: 
        # while g<=GEN: 
            print("-- Generation %i --" % g)
            g += 1
            # Select next generation 
            offspring = toolbox.select(pop, len(pop))
            offspring = toolbox.map(toolbox.clone, offspring)
            
            # Crossover and mutation 
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB: 
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values 
            
            
            # Introduce a rare variant
            if g == 100: 
                
                ind_to_mutate = random.randrange(0,POP)
                toolbox.mutate(offspring[ind_to_mutate])
                mutation_retained = True
                mutation_retained_generations = 0
                del offspring[ind_to_mutate].fitness.values
                
            # Evaluate 
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            pgs = list(toolbox.map(toolbox.evaluate, invalid_ind))
            fitnesses = functions.calc_fitness(pgs, starting_pop_mean, starting_pop_dev)
            starting_fitness_distributions.append(fitnesses)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit,
                
            pop = offspring
                
            heterozygous = np.ndarray(shape=(POP,len(trait_table)))
            for i in range(len(pop)):
                ind = pop[i]
                # Calculate heterozygosity at each locus 
                for x in range(min(len(ind[0]), len(ind[1]))):
                    if (ind[0][x] and not ind[1][x]) or (ind[1][x] and not ind[0][x]): 
                        heterozygous[i][x] = True
                    else: heterozygous[i][x] = False
                    
            
            heterozygous = np.sum(heterozygous,axis=0) / POP
            heterozygosity.append(1 - ((1/len(trait_table))*np.sum(heterozygous)))
            
            # Check if rare mutation retained
            if g >= 100:
                # if any([len(ind[0]) > len(trait_table) for ind in pop]) or any([len(ind[1]) > len(trait_table) for ind in pop]):
                if any([len(ind[0]) > len(trait_table) or len(ind[1]) > len(trait_table) for ind in pop]):
                    "Mutation retained"
                    mutation_retained_generations += 1 
                    # break
                else: 
                    print("Rare mutation lost. Ending simulation")
                    mutation_retained = False
            
                    
            # Some stats
            # if g >= 100: 
                # for ind in pop: 
                #     print(ind.fitness.values)
            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            meanfit = sum(fits) / length
            mean_fits.append(meanfit)
            print("Mean Fit: {0}".format(meanfit), "Max fit: ", min(fits))
            # sum2 = sum(x*x for x in fits)
            # std = abs(sum2 / length - mean**2)**0.5

            # print("  Min %s" % min(fits))
            # print("  Max %s" % max(fits))
            # print("  Avg %s" % mean)
            # print("  Std %s" % std)
            # if g == 100:
            #     fitnesses.sort()
            #     plt.hist(fitnesses)
            #     plt.title("Starting PRS Distribution") 
            #     plt.show()
            #     plt.close()
        
        # Record how long the mutation was retained for 
        
        ax[0].plot(mean_fits)
        ax[1].plot(heterozygosity)
        
        # plt.plot(mean_heights)
        # plt.plot(heterozygosity, linestyle = 'dashed')
        
        retention.append(mutation_retained_generations)
    
    # for fitnesses in starting_fitness_distributions:
    #     sns.distplot(fitnesses,hist = False, kde = True)
    # plt.xlabel('Effect')
    # plt.ylabel('Density')
    # plt.show()
    # plt.close()    
        
    plt.xlabel('Generations')
    plt.ylabel('Fitness')
    plt.show()
    print(retention)
    
    