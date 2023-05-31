import traits_cleaning as traits
import random 
import numpy as np
import statistics
import scipy
from scipy import stats
import CAD_diseaserisk as CAD

import matplotlib.pyplot as plt

random.seed(0)

def make_individual(container, risk_allele_freqs):
    individual = [[],[]]
    for x in range(len(risk_allele_freqs)):
        if random.random() < risk_allele_freqs[x]:
            individual[0].append(True)
        else: individual[0].append(False)
        if random.random() < risk_allele_freqs[x]:
            individual[1].append(True)
        else: individual[1].append(False)
    return container(individual)

def calc_PGS(effect_sizes, individual): 
    score = 0
    known_traits = len(effect_sizes)
    score = sum(effect_sizes[x] for x in range(known_traits) if individual[0][x])
    score += sum(effect_sizes[x] for x in range(known_traits) if individual[1][x])
    # Add rare mutation effects 
    if len(individual[0]) > known_traits:
        for x in range(len(individual[0])-known_traits):
            score += individual[0][known_traits+x]
    if len(individual[1]) > known_traits:
        for x in range(len(individual[1])-known_traits):
            score += individual[1][known_traits+x]
    return score

def rareMutation(ind):
    print("Rare mutation!!") 
    effect = 0.001
    # effect = traits.traits.sample(1).reset_index(drop=True).loc[0]['Effect']
    ind[0].append(effect)
    return ind

# def selRoulette(individuals,k,fit_attr='fitness'): 
#     fits = [getattr(ind,fit_attr).values[0] for ind in individuals]
#     maxFit = np.amax(fits)
#     normalized = [-(x - maxFit)+1 for x in fits]
#     sum_fits = sum(normalized)
#     if sum_fits == 0:
#             fitness_proportions = [1/k for i in normalized]
#     else: fitness_proportions = [i/sum_fits for i in normalized]
    
#     chosen = random.choices(individuals, weights = fitness_proportions, k=k)
#     return chosen   
    
def cxDiploid(ind1, ind2, indpb):
    child1 = [[x for x in ind1[0]],[x for x in ind1[1]]]
    child2 = [[x for x in ind2[0]],[x for x in ind2[1]]]
    if random.randint(0,1) == 1: # Independent assortment
        size = min(len(ind1[0]), len(ind2[0]))
        for i in range(size):
            if random.random() < indpb:
                child1[0][i], child2[0][i] = ind2[0][i], ind1[0][i]
        
        size = min(len(ind1[1]), len(ind2[1]))
        for i in range(size):
            if random.random() < indpb:
                child1[1][i], child2[1][i] = ind2[1][i], ind1[1][i]
    else: 
        size = min(len(ind1[0]), len(ind2[1]))
        for i in range(size):
            if random.random() < indpb:
                child1[0][i], child2[0][i] = ind2[1][i], ind1[1][i]
        
        size = min(len(ind1[1]), len(ind2[0]))
        for i in range(size):
            if random.random() < indpb:
                child1[1][i], child2[1][i] = ind2[0][i], ind1[0][i]
        
    return child1,child2

def selRoulette(individuals, k, fit_attr='fitness'):
    fits = [getattr(ind,fit_attr).values[0] for ind in individuals]
    sum_fits = sum(fits)
    fitness_proportions = [i/sum_fits for i in fits]
    chosen = random.choices(individuals, weights=fitness_proportions, k=k)
    return chosen

def calc_fitness(pgs, mean, dev): 
    cdf = stats.norm.cdf(pgs, loc=mean, scale = dev)
    probability = [stats.norm.cdf(x, mean, dev) for x in pgs]
    risk = [np.interp(x, CAD.percentile, CAD.disease_risk) for x in probability]
    fitness = [1 - x for x in risk]
    return fitness
    
    