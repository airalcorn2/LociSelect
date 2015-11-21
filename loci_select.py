import copy
import csv
import math
import random
import re
import time


def get_individuals(filename):
    """
    Returns a list of all the different individuals contained in the file.
    :param filename:
    :return:
    """
    f = open(filename)
    
    individuals = set()
    
    for line in f:
        
        if "//" in line:
            continue
        
        (individual, sequence) = line[1:].split()
        individuals.add(individual)
        
    return list(individuals)


def get_species_counts(individuals):
    """
    Returns a dictionary with species name as the key and the number of
    individuals in that species as the value.
    :param individuals:
    :return:
    """
    species_counts = {}
    
    for individual in individuals:
        
        species = individual.split("_")[-1]
        if species not in species_counts:
            species_counts[species] = 0
        species_counts[species] += 1
    
    return species_counts


def get_better_loci(filename, cutoff):
    """
    Returns a subset of loci such that each locus includes at least "cutoff"
    different species.
    :param filename:
    :param cutoff:
    :return:
    """
    f = open(filename)
    
    content = f.read()
    f.close()
    loci = re.split(r'//.*|', content)
    
    better_loci = []
    for locus in loci:
        
        found_species = set()
        
        for line in locus.strip().split("\n"):
            
            if line == "":
                continue
           
            (individual, sequence) = line[1:].split()
            found_species.add(individual.split("_")[-1])
        
        if len(found_species) >= cutoff:
            better_loci.append(locus)
    
    return better_loci


def get_loci_score(state, loci_w, data_w, species_w, better_loci,
                   species_counts, total_individuals, total_species,
                   individuals):
    """
    Scoring function with user-specified weights.
    :param state:
    :param loci_w: the included proportion of loci from the original data set (higher is better).
    :param data_w: 1 - the proportion of missing data for the selected loci (higher is better).
    :param species_w: the average proportion of species represented per locus (higher is better).
    :param better_loci:
    :param species_counts:
    :param total_individuals:
    :param total_species:
    :param individuals:
    :return:
    """
    num_loci = sum(state)
    species_loci_counts = {species: 0 for species in species_counts}
    individual_count = 0
    missing_counts = {individual: 0 for individual in individuals}
    total_loci = len(better_loci)
    
    for i in range(total_loci):
        
        if state[i] == 0:
            continue
        
        found_species = set()
        found_individuals = set()
        
        lines = better_loci[i].split("\n")
        for line in lines:
            
            if line == "":
                continue
            
            (individual, sequence) = line[1:].split()
            found_individuals.add(individual)
            individual_count += 1
            
            species = individual.split("_")[-1]
            found_species.add(species)
            
        for species in found_species:
            species_loci_counts[species] += 1
        
        # Keep track of the amount of missing data for each individual.
        for individual in individuals:
            if individual not in found_individuals:
                missing_counts[individual] += 1
    
    num_missing = num_loci * total_individuals - individual_count
    score_comps = [loci_w * float(num_loci) / float(total_loci),
                   data_w * (1 - float(num_missing) / float(num_loci * total_individuals)),
                   species_w * float(sum([species_loci_counts[species] for species in species_loci_counts])) / (float(num_loci) * float(total_species))]
    
    return score_comps, missing_counts


def simulated_annealing(init, temp, cool, iters, loci_w, data_w, species_w,
                        better_loci, species_counts, total_individuals,
                        total_species, individuals):
    """
    Performs the simulated annealing algorithm.
    :param init: specifies the proportion of loci to randomly select for the initial state.
        For example, 0.5 would randomly select approximately 50% of loci for the initial state.
    :param temp: the initial temperature for the algorithm and "cool".
    :param cool: the cooling rate
    :param iters: the number of iterations to run the algorithm.
    :param loci_w:
    :param data_w:
    :param species_w:
    :param better_loci:
    :param species_counts:
    :param total_individuals:
    :param total_species:
    :param individuals:
    :return:
    """
    total_w = loci_w + data_w + species_w
    (loci_w, data_w, species_w) = (loci_w / total_w, data_w / total_w, species_w / total_w)
    proportion = init
    current_state = [1 if random.random() < proportion else 0 for i in range(len(better_loci))]
    best_state = copy.copy(current_state)
    current_score = sum(get_loci_score(current_state, loci_w, data_w, species_w,
                                       better_loci, species_counts,
                                       total_individuals, total_species,
                                       individuals)[0])
    best_score = current_score
    
    iter = 0
    
    while temp > 0.00001 and iter < iters:
        
        if iter % 1000 == 0:
            print("\nIteration: {0}".format(iter))
            print("Temperature: {0:.4f}\n".format(temp))
            print("Current Score: {0:.4f}".format(current_score))
            score_comps, missing_counts = get_loci_score(current_state, loci_w,
                                                         data_w, species_w,
                                                         better_loci, species_counts,
                                                         total_individuals, total_species,
                                                         individuals)
            print("\nProportion of loci used in current state: {0:.4f}".format(score_comps[0] / loci_w))
            print("Proportion of missing data in current state: {0:.4f}".format(1.0 - score_comps[1] / data_w))
            print("Average proportion of species represented per locus in current state: {0:.4f}\n".format(score_comps[2] / species_w))
        
        iter += 1
        temp *= cool
        
        change_index = random.randint(0, len(current_state) - 1)
        candidate_state = copy.copy(current_state)
        candidate_state[change_index] = int(not candidate_state[change_index])
        candidate_score = sum(get_loci_score(candidate_state, loci_w, data_w,
                                             species_w, better_loci, species_counts,
                                             total_individuals, total_species,
                                             individuals)[0])
        
        accept_prob = math.exp((candidate_score - current_score) / temp)
        
        if random.random() < accept_prob:
            
            current_state = candidate_state
            current_score = candidate_score
            
            if current_score > best_score:
                best_score = current_score
                print("{0:.4f}".format(best_score))
                best_state = copy.copy(current_state)
    
    return best_state


def main(test = False):
    # Parameters
    filename = "alpheus_c2.loci"
    cutoff = 4 # Only consider loci with at least "cutoff" species.
    init = 0.5 # Must be between 0 and 1.
    T = 0.1
    cool = 0.999 # Must be between 0 and 1.
    iters = 1e7
    loci_w = 1.0
    data_w = 1.0
    species_w = 1.0
    
    print("Fetching individuals...")
    individuals = get_individuals(filename)
    total_individuals = len(individuals)
    print("{0} different individuals found.\n".format(total_individuals))
    
    print("Counting number of individuals per species...")
    species_counts = get_species_counts(individuals)
    total_species = len(species_counts)
    print("{0} different species found.\n".format(total_species))
    proportions = {species: float(species_counts[species]) / float(total_individuals) for species in species_counts}
    
    print("Selecting loci with at least {0} individuals...".format(cutoff))
    better_loci = get_better_loci(filename, cutoff)
    print("{0} loci found.\n".format(len(better_loci)))
    
    print("Begin simulated annealing...")
    
    best_state = simulated_annealing(init, T, cool, iters, loci_w, data_w,
                                     species_w, better_loci, species_counts, total_individuals,
                                     total_species, individuals)
    print("\nSimulated annealing complete.")
    
    score_comps, missing_counts = get_loci_score(best_state, loci_w, data_w, species_w,
                                                 better_loci, species_counts,
                                                 total_individuals, total_species,
                                                 individuals)
    print("\nProportion of loci used in best discovered state: {0:.4f}".format(score_comps[0] / loci_w))
    print("Proportion of missing data in best discovered state: {0:.4f}".format(1.0 - score_comps[1] / data_w))
    print("Average proportion of species represented per locus in current state: {0:.4f}".format(score_comps[2] / species_w))
    
    output = open("selected_loci.txt", "w")
    
    for i in range(len(better_loci)):
        if better_loci[i]:
            output.write(str(i) + "\n")

    if test:
        run_test(T, cool, iters, better_loci, species_counts, total_individuals,
                 total_species, individuals)


def run_test(T, cool, iters, better_loci, species_counts, total_individuals,
             total_species, individuals):
    """
    Test the program.
    :param T:
    :param cool:
    :param iters:
    :param better_loci:
    :param species_counts:
    :param total_individuals:
    :param total_species:
    :param individuals:
    :return:
    """
    reader = csv.DictReader(open("parameters.csv"))
    output = open("results.csv", "w")
    fieldnames = reader.fieldnames + ["score", "loci", "data", "species"]
    writer = csv.DictWriter(output, fieldnames)
    writer.writeheader()
    
    start = time.time()
    
    for row in reader:
        init = float(row["init"])
        loci_w = float(row["loci_w"])
        data_w = float(row["data_w"])
        species_w = float(row["species_w"])
        best_state = simulated_annealing(init, T, cool, iters, loci_w, data_w, species_w)
        score_comps, missing_counts = get_loci_score(best_state, loci_w, data_w, species_w,
                                                     better_loci, species_counts,
                                                     total_individuals, total_species,
                                                     individuals)
        new_row = {col: row[col] for col in row}
        new_row["score"] = sum(score_comps)
        new_row["loci"] = score_comps[0] / loci_w
        new_row["data"] = score_comps[1] / data_w
        new_row["species"] = score_comps[2] / species_w
        writer.writerow(new_row)
    
    output.close()
    
    end = time.time()
    
    print("Total time: {0} s".format(end - start))


if __name__ == "__main__":
    main()