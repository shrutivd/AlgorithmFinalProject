import time
import matplotlib.pyplot as plt
from itertools import permutations, combinations
import random
import numpy as np

x = [0,30,60,70,150,100,160,50,80,15]
y = [1,20,10,45,-10,25,110,60,90,120]
genes = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
#,"OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV", "XX", "YY", , "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM","NN"
# , "U", "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"
city_count = len(genes)
print("city count :", city_count)
#x = random.sample(range(0, 150), city_count)
#print("x values :", x)
#y = random.sample(range(0, 10), city_count)
#print("y values :", y)
geneCoordinates = dict(zip(genes, zip(x, y)))



# fig, ax = plt.subplots()

# ax.grid(False)  # Grid

# for i, (city, (city_x, city_y)) in enumerate(city_coords.items()):
#     icon = city
#     ax.scatter(city_x, city_y, s=1200, zorder=2)
#     ax.annotate(icon, (city_x, city_y), fontsize=12, ha='center', va='center', zorder=3)
#
#     # Connect cities with opaque lines
#     for j, (other_city, (other_x, other_y)) in enumerate(city_coords.items()):
#         if i != j:
#             ax.plot([city_x, other_x], [city_y, other_y], color='gray', linestyle='-', linewidth=1, alpha=0.1)
#
# fig.set_size_inches(16, 12)
# plt.show()


def generatePopulation(cities_list, n_population=250):
    population = []

    for i in range(n_population):
        chromosome = genes.copy()
        random.shuffle(chromosome)
        population.append(chromosome.copy())

    return population


def distanceBetweenGenes(gene1Name, gene2Name):
    gene1Coordinates = geneCoordinates[gene1Name]
    gene2Coordinates = geneCoordinates[gene2Name]
    return np.sqrt(np.sum((np.array(gene1Coordinates) - np.array(gene2Coordinates)) ** 2))


def totalChromosomeDistance(chromosome):
    chromosomeDistance = 0
    for i in range(0, len(chromosome)):
        if (i == len(chromosome) - 1):
            chromosomeDistance += distanceBetweenGenes(chromosome[i], chromosome[0])
        else:
            chromosomeDistance += distanceBetweenGenes(chromosome[i], chromosome[i + 1])
    return chromosomeDistance


def fitnessProbability(population):
    totalPopulationDistance = []
    for i in range(0, len(population)):
        totalPopulationDistance.append(totalChromosomeDistance(population[i]))

    maxCost = max(totalPopulationDistance)
    populationFitness = maxCost - totalPopulationDistance
    populationFitnessSum = sum(populationFitness)
    populationFitnessProbability = populationFitness / populationFitnessSum
    return populationFitnessProbability


def rouletteWheel(population, fitnessProbability):
    # population fitness probability = PFP
    PFP_cummulativeSum = fitnessProbability.cumsum()
    booleanProbabilityArray = PFP_cummulativeSum < np.random.uniform(0, 1, 1)
    selectedChoromosomeIndex = len(booleanProbabilityArray[booleanProbabilityArray == True]) - 1
    return population[selectedChoromosomeIndex]


def crossover(parent1, parent2):
    n_cities_cut = len(genes) - 1
    cut = round(random.uniform(1, n_cities_cut))
    offspring1 = []
    offspring2 = []

    offspring1 = parent1[0:cut]
    offspring1 += [gene for gene in parent2 if gene not in offspring1]

    offspring2 = parent2[0:cut]
    offspring2 += [gene for gene in parent1 if gene not in offspring2]

    return offspring1, offspring2


def mutation(offspring):
    n_cities_cut = len(genes) - 1
    index1 = round(random.uniform(0, n_cities_cut))
    index2 = round(random.uniform(0, n_cities_cut))

    temp = offspring[index1]
    offspring[index1] = offspring[index2]
    offspring[index2] = temp
    return(offspring)


def run_ga(geneName, n_population, cut_percentage, crossoverPercentage, mutationPercentage):
    population = generatePopulation(geneName, n_population)
    fitnessProba = fitnessProbability(population)

    parentsList = []
    for i in range(0, int(crossoverPercentage * n_population)):
        parentsList.append(rouletteWheel(population, fitnessProba))

    offspringList = []
    for i in range(0, len(parentsList), 2):
        offspring1, offspring2 = crossover(parentsList[i], parentsList[i + 1])

        mutateThreashold = random.random()
        if (mutateThreashold > (1 - mutationPercentage)):
            offspring1 = mutation(offspring1)
        #         print("Offspring 1 mutated", offspring_1)

        mutateThreashold = random.random()
        if (mutateThreashold > (1 - mutationPercentage)):
            offspring2 = mutation(offspring2)
        #         print("Offspring 2 mutated", offspring_2)

        offspringList.append(offspring1)
        offspringList.append(offspring2)

    mixedOffspring = parentsList + offspringList

    fitnessProba = fitnessProbability(mixedOffspring)  # Corrected variable name
    sortedFitnessIndices = np.argsort(fitnessProba)[::-1]
    bestFitnessIndices = sortedFitnessIndices[0:n_population]
    bestMixedOffsrping = []
    for i in bestFitnessIndices:
        bestMixedOffsrping.append(mixedOffspring[i])

    prevMinimumDistance = totalChromosomeDistance(bestMixedOffsrping[0])
    curMinimumDistance = prevMinimumDistance * 1000
    # print(prev_minimum_distance, cur_minimum_distance)
    generations = [1]
    distances = [prevMinimumDistance]
    gen = 1
    index_minimum = 0
    while (abs(curMinimumDistance-prevMinimumDistance) > cut_percentage * curMinimumDistance):
        # if (i%10 == 0):
        prevMinimumDistance = curMinimumDistance
        distance = 0
        for _ in range(20):
            # print("Generation: ", generations)
            gen += 1

            fitnessProba = fitnessProbability(bestMixedOffsrping)  # Corrected variable name
            parentsList = []
            for i in range(0, int(crossoverPercentage * n_population)):
                parentsList.append(rouletteWheel(bestMixedOffsrping, fitnessProba))

            offspringList = []
            for i in range(0, len(parentsList), 2):
                offspring1, offspring_2 = crossover(parentsList[i], parentsList[i + 1])

                mutateThreashold = random.random()
                if (mutateThreashold > (1 - mutationPercentage)):
                    offspring1 = mutation(offspring1)

                mutateThreashold = random.random()
                if (mutateThreashold > (1 - mutationPercentage)):
                    offspring_2 = mutation(offspring_2)

                offspringList.append(offspring1)
                offspringList.append(offspring_2)

            mixedOffspring = parentsList + offspringList
            fitnessProba = fitnessProbability(mixedOffspring)  # Corrected variable name
            sortedFitnessIndices = np.argsort(fitnessProba)[::-1]
            bestFitnessIndices = sortedFitnessIndices[0:int(0.8 * n_population)]

            bestMixedOffsrping = []
            for i in bestFitnessIndices:
                bestMixedOffsrping.append(mixedOffspring[i])

            oldPopulationIndices = [random.randint(0, (n_population - 1)) for j in range(int(0.2 * n_population))]
            for i in oldPopulationIndices:
                #             print(i)
                bestMixedOffsrping.append(population[i])

            random.shuffle(bestMixedOffsrping)

            totalChromosomeDistances = []
            for i in range(0, n_population):
                totalChromosomeDistances.append(totalChromosomeDistance(bestMixedOffsrping[i]))

            index_minimum = np.argmin(totalChromosomeDistances)
            distance += min(totalChromosomeDistances)

        curMinimumDistance = distance/20
        generations.append(gen)
        distances.append(curMinimumDistance)

    return bestMixedOffsrping, generations, index_minimum, distances


n_population = 50
crossoverPercentage = 0.8
mutationPercentage = 0.2
cut_threshold = 0.001
distances = []
averageGenerations = 0
averageDistance = 0
for _ in range(5):
    bestMixedOffsrping, generations, index_minimum, distances = run_ga(genes, n_population, cut_threshold, crossoverPercentage, mutationPercentage)
    averageDistance += distances[-1]
    averageGenerations += generations[-1]

averageDistance = averageDistance / 5
averageGenerations = averageGenerations / 5
print("Average Generations: ", averageGenerations, ", Minimum distance: ", averageDistance)
plt.plot(generations, distances)
plt.show()

shortest_path = bestMixedOffsrping[index_minimum]

x_shortest = []
y_shortest = []
for city in shortest_path:
    x_value, y_value = geneCoordinates[city]
    x_shortest.append(x_value)
    y_shortest.append(y_value)

x_shortest.append(x_shortest[0])
y_shortest.append(y_shortest[0])

fig, ax = plt.subplots()
ax.plot(x_shortest, y_shortest, '--go', label='Best Route', linewidth=2.5)
plt.legend()

for i in range(len(x)):
    for j in range(i + 1, len(x)):
        ax.plot([x[i], x[j]], [y[i], y[j]], 'k-', alpha=0.09, linewidth=1)

plt.title(label="TSP Best Route Using GA",
          fontsize=25,
          color="k")

str_params = '\n' + str(generations) + ' Generations\n' + str(n_population) + ' Population Size\n' + str(
    crossoverPercentage) + ' Crossover\n' + str(mutationPercentage) + ' Mutation'
plt.suptitle("Total Distance Travelled: " +
             str(round(distances[-1], 3)) +
             str_params, fontsize=18, y=1.047)

for i, txt in enumerate(shortest_path):
    ax.annotate(str(i + 1) + "- " + txt, (x_shortest[i], y_shortest[i]), fontsize=20)

fig.set_size_inches(16, 12)
# plt.grid(color='k', linestyle='dotted')
plt.savefig('solution.png')
plt.show()