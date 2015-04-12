import copy
import csv
import math
import random
import re
import time


# Returns a list of all the different individuals contained in the file.
def getIndividuals(filename):
    
    f = open(filename)
    
    individuals = set()
    
    for line in f:
        
        if "//" in line:
            continue
        
        (individual, sequence) = line[1:].split()
        individuals.add(individual)
        
    return list(individuals)


# Returns a dictionary with species name as the key and the number of individuals in that species as the value.
def getSpeciesCounts(individuals):
    
    speciesCounts = {}
    
    for individual in individuals:
        
        species = individual.split("_")[-1]
        if species not in speciesCounts:
            speciesCounts[species] = 0
        speciesCounts[species] += 1
    
    return speciesCounts


# Returns a subset of loci such that each locus includes at least "cutoff" different species.
def getBetterLoci(filename, cutoff):
    
    f = open(filename)
    
    content = f.read()
    f.close()
    loci = re.split(r'//.*|', content)
    
    betterLoci = []
    for locus in loci:
        
        foundSpecies = set()
        
        for line in locus.strip().split("\n"):
            
            if line == "":
                continue
           
            (individual, sequence) = line[1:].split()
            foundSpecies.add(individual.split("_")[-1])
        
        if len(foundSpecies) >= cutoff:
            betterLoci.append(locus)
    
    return betterLoci


# Scoring function with user-specified weights. Includes three components:
# 1. The included proportion of loci from the original data set (higher is better). (weight = "lociW")
# 2. 1 - the proportion of missing data for the selected loci (higher is better). (weight = "dataW")
# 3. The average proportion of species represented per locus (higher is better). (weight = "speciesW")
def getLociScore(state, lociW, dataW, speciesW):
    
    global betterLoci, speciesCounts, totalIndividuals, totalSpecies, individuals
    
    numLoci = sum(state)
    speciesLociCounts = {species: 0 for species in speciesCounts.keys()}
    individualCount = 0
    missingCounts = {individual: 0 for individual in individuals}
    totalLoci = len(betterLoci)
    
    for i in range(totalLoci):
        
        if state[i] == 0:
            continue
        
        foundSpecies = set()
        foundIndividauals = set()
        
        lines = betterLoci[i].split("\n")
        for line in lines:
            
            if line == "":
                continue
            
            (individual, sequence) = line[1:].split()
            foundIndividauals.add(individual)
            individualCount += 1
            
            species = individual.split("_")[-1]
            foundSpecies.add(species)
            
        for species in foundSpecies:
            speciesLociCounts[species] += 1
        
        # Keep track of the amount of missing data for each individual.
        for individual in individuals:
            if individual not in foundIndividauals:
                missingCounts[individual] += 1
    
    numMissing = numLoci * totalIndividuals - individualCount
    scoreComps = [lociW * float(numLoci) / float(totalLoci),
                  dataW * (1 - float(numMissing) / float(numLoci * totalIndividuals)),
                  speciesW * float(sum([speciesLociCounts[species] for species in speciesLociCounts.keys()])) / (float(numLoci) * float(totalSpecies))]
    
    return scoreComps, missingCounts


# Performs simulated annealing algorithm.
# "init" specifies the proportion of loci to randomly select for the initial state.
#       For example, 0.5 would randomly select approximately 50% of loci for the initial state.
# "temp" is the initial temperature for the algorithm and "cool" is the cooling rate.
# "iters" is the number of iterations to run the algorithm.
def simulatedAnnealing(init, temp, cool, iters, lociW, dataW, speciesW):
    
    global betterLoci
    
    totalW = lociW + dataW + speciesW
    (lociW, dataW, speciesW) = (lociW / totalW, dataW / totalW, speciesW / totalW)
    proportion = init
    currentState = [1 if random.random() < proportion else 0 for i in range(len(betterLoci))]
    bestState = copy.copy(currentState)
    currentScore = sum(getLociScore(currentState, lociW, dataW, speciesW)[0])
    bestScore = currentScore
    
    iter = 0
    
    while temp > 0.00001 and iter < iters:
        
        if iter % 1000 == 0:
            print("\nIteration: {0}".format(iter))
            print("Temperature: {0}\n".format(temp))
            print("Current Score: {0}".format(currentScore))
            scoreComps, missingCounts = getLociScore(currentState, lociW, dataW, speciesW)
            print("\nProportion of loci used in current state: {0}".format(scoreComps[0] / lociW))
            print("Proportion of missing data in current state: {0}".format(1.0 - scoreComps[1] / dataW))
            print("Average proportion of species represented per locus in current state: {0}".format(scoreComps[2] / speciesW))
            print("")
        
        iter += 1
        temp *= cool
        
        changeIndex = random.randint(0, len(currentState) - 1)
        candidateState = copy.copy(currentState)
        candidateState[changeIndex] = int(not candidateState[changeIndex])
        candidateScore = sum(getLociScore(candidateState, lociW, dataW, speciesW)[0])
        
        acceptProb = math.exp((candidateScore - currentScore) / temp)
        
        if random.random() < acceptProb:
            
            currentState = candidateState
            currentScore = candidateScore
            
            if currentScore > bestScore:
                bestScore = currentScore
                print(bestScore)
                bestState = copy.copy(currentState)
    
    return bestState


# Parameters
filename = "alpheus_c2.loci"
cutoff = 4 # Only consider loci with at least "cutoff" species.
init = 0.5 # Must be between 0 and 1.
T = 0.1
cool = 0.999 # Must be between 0 and 1.
iters = 1e7
lociW = 1.0
dataW = 1.0
speciesW = 1.0

print("Fetching individuals...")
individuals = getIndividuals(filename)
totalIndividuals = len(individuals)
print("{0} different individuals found.\n".format(totalIndividuals))

print("Counting number of individuals per species...")
speciesCounts = getSpeciesCounts(individuals)
totalSpecies = len(speciesCounts.keys())
print("{0} different species found.\n".format(totalSpecies))
proportions = {species: float(speciesCounts[species]) / float(totalIndividuals) for species in speciesCounts.keys()}

print("Selecting loci with at least {0} individuals...".format(cutoff))
betterLoci = getBetterLoci(filename, cutoff)
print("{0} loci found.\n".format(len(betterLoci)))

print("Begin simulated annealing...")

bestState = simulatedAnnealing(init, T, cool, iters, lociW, dataW, speciesW)
print("\nSimulated annealing complete.")

scoreComps, missingCounts = getLociScore(bestState, lociW, dataW, speciesW)
print("\nProportion of loci used in best discovered state: {0}".format(scoreComps[0] / lociW))
print("Proportion of missing data in best discovered state: {0}".format(1.0 - scoreComps[1] / dataW))
print("Average proportion of species represented per locus in current state: {0}".format(scoreComps[2] / speciesW))

'''
# Testing
reader = csv.DictReader(open("parameters.csv"))
output = open("results.csv", "w")
fieldnames = reader.fieldnames + ["score", "loci", "data", "species"]
writer = csv.DictWriter(output, fieldnames)
writer.writeheader()

start = time.time()

for row in reader:
    init = float(row["init"])
    lociW = float(row["lociW"])
    dataW = float(row["dataW"])
    speciesW = float(row["speciesW"])
    bestState = simulatedAnnealing(init, T, cool, iters, lociW, dataW, speciesW)
    scoreComps, missingCounts = getLociScore(bestState, lociW, dataW, speciesW)
    newRow = {col: row[col] for col in row.keys()}
    newRow["score"] = sum(scoreComps)
    newRow["loci"] = scoreComps[0] / lociW
    newRow["data"] = scoreComps[1] / dataW
    newRow["species"] = scoreComps[2] / speciesW
    writer.writerow(newRow)

output.close()

end = time.time()

print("Total time: {0} s".format(end - start))
'''
