#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sstream>
#include <algorithm>
#include <array>
#include <iomanip>
#include <bitset>
#include <stack>
#include <cmath>
template<size_t bit_len>
class Solution {
    public:
        Solution(const Solution& sol): i(sol.i), length(sol.length), value(sol.value), weight(sol.weight) {
            this->chromosome = sol.chromosome;
        }
        Solution& operator=(const Solution& sol) {
            if (this == &sol) {
                return *this;
            }
            this->i = sol.i;
            this->length = sol.length;
            this->value = sol.value;
            this->weight = sol.weight;
            this->chromosome = sol.chromosome;
            return *this;
        }
        Solution(){};
        ~Solution() {}
        void setChromosome(std::bitset<bit_len> chromosome) {
            this->chromosome = chromosome;
            this->value = 0;
            this->weight = 0;
        }
        void print() {
            for (int i = 0; i < length; i++) {
                if (chromosome[i]) {
                    std::cout << i << " ";
                }
            }
            std::cout << " [V: " << value << "]";
            std::cout << "[W: " << weight << "]\n";
        }
        int i;
        // int *chromosome = nullptr;
        std::bitset<bit_len> chromosome;
        int length = bit_len;
        double value;
        double weight;
};
// template<size_t bit_len>

const size_t TODO = 111;
class GA {
    public:
        std::string filename;
        bool debug = false;
        unsigned int seed = 0;
        GA(
            std::vector<std::array<double, 3>> items, 
            int length, 
            int knapsackMax, 
            int populationSize, 
            int generations, 
            double crossoverRate, 
            double mutationRate, 
            int elitismCount,
            bool memetic,
            bool debug = false,
            std::string filename = ""
        ) : filename(filename),
            numItems(length), 
            knapsackMax(knapsackMax), 
            populationSize(populationSize), 
            generations(generations), 
            crossoverRate(crossoverRate), 
            mutationRate(mutationRate), 
            elitismCount(elitismCount),
            memetic(memetic),
            debug(debug)
        {
            if (debug) {
                std::cout << "Debugging" << std::endl;
            }
            for (char c : filename) {
                seed += c;
            }
            // std::cout << "Seed value for " << filename << ": " << seed << std::endl;
            if (seed > 0) {
                srand(seed);
            }
            // std::cout << "GA constructor" << std::endl;
            // std::cout << "Population created" << std::endl;
            this->items = new std::array<double, 3>[length];
            for (int i = 0; i < length; i++) {
                this->items[i] = items[i];
                this->avgWeight += items[i][1];
                if (debug) 
                    std::cout << "Item " << i << ": " << items[i][0] << " " << items[i][1] << " " << items[i][2] << std::endl;
                
            }
            this->avgWeight /= length;

            this->itemsEfficiencyIndex = new int[length];
            for (int i = 0; i < length; i++) {
                this->itemsEfficiencyIndex[i] = i;
            }
            std::sort(this->itemsEfficiencyIndex, this->itemsEfficiencyIndex + length, 
                [this](int a, int b) {
                    return this->items[a][2] > this->items[b][2];
                }
            );
            for (int i = 0; i < length; i++) {
                if (debug) 
                    std::cout << "Item " << i << ": " << items[this->itemsEfficiencyIndex[i]][0] << " " << items[this->itemsEfficiencyIndex[i]][1] << " " << items[this->itemsEfficiencyIndex[i]][2] << std::endl;
                
            }
            // std::cout << "Items populated" << std::endl;
            this->init();
            // std::cout << "Init complete" << std::endl;
            this->bestFitnessPerGeneration = new double[generations];

            int solSize = round(knapsackMax / avgWeight);
            long long combinations = numItems;
            for (int i = 1; i < solSize; i++) {
                combinations *= numItems - i;
                combinations /= i + 1;
            }
            double ratio = (double)combinations / populationSize;
            double logRatio = log10(ratio);

            this->mutationRate = 0.01 + (0.1 - 0.01) * (logRatio - 0) / (log10(combinations) - 0);
        }
        ~GA(){
            // std::cout << "GA destructor" << std::endl;
            delete [] items;
            delete [] itemsEfficiencyIndex;
            delete [] bestFitnessPerGeneration;
        }
        void run(){
            if (debug) {
                std::cout << "Running GA" << std::endl;
            }
            double bestFitness = 0;
            for (int i = 0; i < generations; i++) {
                if (debug)
                    std::cout << "Generation " << i << std::endl;
                eval();
                if (population[0].value > bestFitness) {
                    if (debug) {
                        std::cout << "Generation " << i << " new best: ";
                        population[0].print();
                    }
                    bestFitness = population[0].value;
                    generationOfBest = i;
                }
                if (debug)
                    std::cout << "Evaluated" << std::endl;
                // for (int i = 0; i < populationSize; i++) {
                //     std::cout << "\t" << i << ": ";
                //     population[i].print();
                // }
                // std::cout << "Total fitness: " << totalFitness << std::endl;
                selectParents();
                if (debug)
                    std::cout << "Parents selected: " << std::endl;
                // for (int i = 0; i < parents.size(); i++) {
                //     std::cout << "\t" << i << ": ";
                //     population[parents[i]].print();
                // }
                // std::cout << std::endl;
                crossover();
                if (debug)
                    std::cout << "Crossover complete" << std::endl;
                // std::cout << "Offspring: " << std::endl;
                // for (int i = 0; i < offspring.size(); i++) {
                //     std::cout << "\t" << i << ": ";
                //     offspring[i].print();
                // }
                mutate();
                if (debug)
                    std::cout << "Mutation complete" << std::endl;
                // std::cout << "Offspring: " << std::endl;
                // for (int i = 0; i < offspring.size(); i++) {
                //     std::cout << "\t" << i << ": ";
                //     offspring[i].print();
                // }
                if (memetic) {
                    for (int i = 0; i < offspring.size(); i++) {
                        localSearch(offspring[i]);
                    }
                }
                replace();
                bestFitnessPerGeneration[i] = population[0].value;
                if (debug)
                    std::cout << "Replace complete" << std::endl;

            }
            // print();
            // std::cout << "Final population: " << std::endl;
            // for (int i = 0; i < populationSize; i++) {
            //     std::cout << "\t" << i << ": ";
            //     population[i].print();
            // }
            // std::cout << "\t Gen " << std::setw(4) << std::setfill('0') << generationOfBest << ". Best sol: ";
            // population[0].print();
            // std::cout << "Value: " << fitness(population[0].chromosome) << std::endl;
            // std::cout << "Done" << std::endl;
        }
        std::vector<Solution<111>> population;
        std::vector<Solution<111>> offspring;
        int knapsackMax;
        int populationSize;
        int generations;
        double crossoverRate;
        double mutationRate;
        int elitismCount;
        bool memetic;
        std::vector<int> parents;
        // pair->first is id, pair->second is weight
        std::array<double, 3> *items;
        int numItems;
        double avgWeight = 0;
        double totalFitness;
        int generationOfBest;
        int* itemsEfficiencyIndex;
        double *bestFitnessPerGeneration;
    private:
        // pair->first is chromosome, pair->second is fitness
        // Solution *population;
        void init() {
            if (debug) {
                std::cout << "Initialising population" << std::endl;
                std::cout << "numItems: " << numItems << std::endl;
                std::cout << "knapsackMax: " << knapsackMax << std::endl;
            }
            int consecutiveDuplicates = 0;
            for (int i = 0; i < populationSize && consecutiveDuplicates < 1000; i++) {
                std::bitset<111> chromosome;
                // std::cout << "Chromosome " << i << std::endl;
                // std::cout << "numItems: " << numItems << std::endl;
                // std::cout << "knapsackMax: " << knapsackMax << std::endl;
                // random chromosome length grater than 2
                // std::cout << "Generating chromosome of length " << chromosomeLength << std::endl;
                int index;
                int counter = 0;
                while (weight(chromosome) < knapsackMax) {
                    // std::cout << "Here 1\n";
                    do {
                        index = rand() % numItems;
                    } while (chromosome[index]);
                    // std::cout << "Here 2\n";
                    ++ counter;
                    // std::cout << "Here 3\n";
                    chromosome.set(index, 1);
                    // if (counter > 1000) {
                    //     std::cout << chromosome << std::endl;
                    //     std:: cout << "Weight: " << weight(chromosome) << std::endl;
                    // }
                }
                // undo last flip that set it over weight limit
                chromosome.flip(index);
                // std::cout << "Here 4\n";
                double weight = this->weight(chromosome);

                if (weight == 0) {
                    // std::cout << "Chromosome " << i << " has weight 0" << std::endl;
                    i--;
                    continue;
                }
                // std::cout << "Here 7\n";
                bool chromosomeExists = exists(chromosome);
                if (chromosomeExists) {
                    // std::cout << "Chromosome " << i << " already exists" << std::endl;
                    ++consecutiveDuplicates;
                    i--;
                    continue;
                }
                Solution<111> sol;
                sol.setChromosome(chromosome);
                sol.weight = weight;
                sol.value = fitness(chromosome);
                sol.length = numItems;
                population.push_back(sol);
                consecutiveDuplicates = 0;
            }
            if (consecutiveDuplicates >= 1000) {
                if (debug)
                    std::cout << "Unable to generate unique population, reducing populationSize to " << population.size() << std::endl;
                populationSize = population.size();
            }
            
            std::sort(population.begin(), population.end(), 
                [](const Solution<111>& a, const Solution<111>& b) {
                    return a.value > b.value;
                }
            );
        }
        bool exists(std::bitset<111> chromosome, std::vector<Solution<111>> *pop = nullptr) {
            if (pop == nullptr) {
                pop = &population;
            }
            for (int i = 0; i < pop->size() && (*pop)[i].length != 0; i++) {
                if (chromosome == (*pop)[i].chromosome) {
                    return true;
                }
                // population is sorted by value
                if (fitness(chromosome) > (*pop)[i].value) {
                    return false;
                }
            }
            return false;
        }
        double fitness(std::bitset<111> chromosome) {
            double value = 0;
            for (int i = 0; i < numItems; i++) {
                if (chromosome[i]) {
                    value += items[i][0];
                }
            }
            return value;
        }
        double weight(std::bitset<111> chromosome) {
            double weight = 0;
            for (int i = 0; i < numItems; i++) {
                if (chromosome[i]) {
                    weight += items[i][1];
                }
            }
            return weight;
        }
        double weight(std::vector<int> items) {
            double weight = 0;
            for (int i = 0; i < items.size(); i++) {
                weight += this->items[items[i]][1];
            }
            return weight;
        }
        double fitness(std::vector<int> items) {
            double value = 0;
            for (int i = 0; i < items.size(); i++) {
                value += this->items[items[i]][0];
            }
            return value;
        }
        void eval(){
            totalFitness = 0;
            for (int i = 0; i < populationSize; i++) {
                totalFitness += population[i].value;
            }
            std::sort(population.begin(), population.end(), 
                [](const Solution<111>& a, const Solution<111>& b) {
                    return a.value > b.value;
                }
            );
        }
        void selectParents(){
            bool selected[populationSize];
            for (int i = 0; i < populationSize; i++) {
                selected[i] = false;
            }
            parents.clear();
            double numParents = populationSize * 0.1;
            if (numParents < 2) {
                numParents = populationSize;
            }
            int numSelected = 0;
            while (numSelected < numParents) {
                int index = rand() % (populationSize - numSelected);
                int counter = 0;
                while (index > 0 && counter < populationSize) {
                    if (selected[counter] == false) {
                        --index;
                    }
                    ++counter;
                }
                while(counter < populationSize && selected[counter] == true) {
                    ++counter;
                }
                if (counter == populationSize) {
                    continue;
                }
                ++numSelected;
                selected[counter] = true;
                parents.push_back(counter);
            }
            // std::cout << "Parents selected: " << std::endl;
            // for (int i = 0; i < parents.size(); i++) {
            //     std::cout << "\tP " << i << ": ";
            //     population[parents[i]].print();
            // }
        }
        void crossover(){
            // std::cout << "Crossover" << std::endl;
            offspring.clear();
            while (parents.size() > 1) {
                int currParents[2];
                currParents[0] = parents.back();
                parents.pop_back();
                currParents[1] = parents.back();
                parents.pop_back();
                // place shortest parent first
                if (population[currParents[0]].length > population[currParents[1]].length) {
                    int temp = currParents[0];
                    currParents[0] = currParents[1];
                    currParents[1] = temp;
                }
                // std::cout << "Parent 1: ";
                // population[currParents[0]].print();
                // std::cout << "Parent 2: ";
                // population[currParents[1]].print();
                double r = (double)rand() / RAND_MAX;
                if (r < crossoverRate) {
                    std::bitset<111> child1;
                    std::bitset<111> child2;
                    std::vector<int> child1Items, child2Items;
                    std::bitset<111> mask;
                    // std::cout << "Mask: " << mask << std::endl;
                    for (int i = 0; i < numItems; i++) {
                        int swap = rand() % 2;
                        mask.set(i, swap);
                        if (swap) {
                            child1.set(i, population[currParents[0]].chromosome[i]);
                            child2.set(i, population[currParents[1]].chromosome[i]);
                        } else {
                            child1.set(i, population[currParents[1]].chromosome[i]);
                            child2.set(i, population[currParents[0]].chromosome[i]);
                        }
                        if (child1[i]) {
                            child1Items.push_back(i);
                        }
                        if (child2[i]) {
                            child2Items.push_back(i);
                        }
                    }
                    if (debug) {
                        std::cout << "Mask: ";
                        for (int i = 0; i < numItems; i++) {
                            if (mask[i])
                                std::cout << "\033[32m" << mask[i] << "\033[0m";
                            else 
                                std::cout << "\033[33m" << mask[i] << "\033[0m";
                        }
                        std::cout << std::endl;
                        std::cout << "P1  : ";
                        for (int i = 0; i < numItems; i++) {
                            std::cout << "\033[32m" << population[currParents[0]].chromosome[i] << "\033[0m";
                        }
                        std::cout << std::endl;
                        std::cout << "P2  : ";
                        for (int i = 0; i < numItems; i++) {
                            std::cout << "\033[33m" << population[currParents[1]].chromosome[i] << "\033[0m";
                        }
                        std::cout << std::endl;
                        std::cout << "C1  : ";
                        for (int i = 0; i < numItems; i++) {
                            if (mask[i])
                                std::cout << "\033[32m";
                            else 
                                std::cout << "\033[33m";
                            std::cout << child1[i] << "\033[0m";
                        }
                        std::cout << std::endl;
                        std::cout << "C2  : ";
                        for (int i = 0; i < numItems; i++) {
                            if (!mask[i])
                                std::cout << "\033[32m";
                            else 
                                std::cout << "\033[33m";
                            std::cout << child2[i] << "\033[0m";
                        }
                    }
                    
                    while (weight(child1) > knapsackMax) {
                        // unset random bit until weight is below knapsackMax
                        int index = rand() % child1Items.size();
                        if (debug) {
                            std::cout << "Unsetting bit " << child1Items[index] << " of child1 to reduce weight" << std::endl;
                            for (int i = 0; i < numItems; i++) {
                                if (i == child1Items[index])
                                    std::cout << "\033[33m";
                                std::cout << child1[i] << "\033[0m";
                            }
                            std::cout << std::endl;
                        }
                        child1.set(child1Items[index], 0);
                        if (debug) {
                            for (int i = 0; i < numItems; i++) {
                                if (i == child1Items[index])
                                    std::cout << "\033[33m";
                                std::cout << child1[i] << "\033[0m";
                            }
                            std::cout << std::endl;
                        }
                        child1Items.erase(child1Items.begin() + index);
                    }
                    while (weight(child2) > knapsackMax) {
                        // unset random bit until weight is below knapsackMax
                        int index = rand() % child2Items.size();
                        if (debug) {
                            std::cout << "Unsetting bit " << child2Items[index] << " of child2 to reduce weight" << std::endl;
                            for (int i = 0; i < numItems; i++) {
                                if (i == child2Items[index])
                                    std::cout << "\033[33m";
                                std::cout << child2[i] << "\033[0m";
                            }
                            std::cout << std::endl;
                        }
                        child2.set(child2Items[index], 0);
                        if (debug) {
                            for (int i = 0; i < numItems; i++) {
                                if (i == child2Items[index])
                                    std::cout << "\033[33m";
                                std::cout << child2[i] << "\033[0m";
                            }
                            std::cout << std::endl;
                        }
                        child2Items.erase(child2Items.begin() + index);
                    }
                    // std::cout << std::endl;
                    Solution<111> sol1;
                    sol1.setChromosome(child1);
                    sol1.value = fitness(child1);
                    sol1.weight = weight(child1);
                    sol1.length = numItems;
                    offspring.push_back(sol1);

                    Solution<111> sol2;
                    sol2.setChromosome(child2);
                    sol2.value = fitness(child2);
                    sol2.weight = weight(child2);
                    sol2.length = numItems;
                    offspring.push_back(sol2);

                } 
                else 
                {
                    // as per Sastry, K., Goldberg, D., & Kendall, G. (2014). Genetic algorithms. In Search methodologies (pp. 96). Springer, Boston, MA.
                    Solution<111> sol1;
                    Solution<111> sol2;
                    sol1 = population[currParents[0]];
                    sol2 = population[currParents[1]];
                    offspring.push_back(sol1);
                    offspring.push_back(sol2);
                }
            }
            if (debug){
                std::cout << "Offspring: " << std::endl;
                for (int i = 0; i < offspring.size(); i++) {
                    std::cout << "\tO " << i << ": ";
                    for (int j = 0; j < offspring[i].length; j++) {
                        if (offspring[i].chromosome[j]) {
                            std::cout << j << " ";
                        }
                    }
                    std::cout << ", Value: " << offspring[i].value;
                    std::cout << ", Weight: " << offspring[i].weight << std::endl;
                }
            }
        }
        void mutate(){
            bool selected[numItems];
            int numTested = 0;
            int index;
            double r;
            for (int i = 0; i < offspring.size(); i++) {
                // if (debug){
                //     std::cout << "Offspring " << i << std::endl;
                //     offspring[i].print();
                // }
                r = (double)rand() / RAND_MAX;
                if (debug) {
                    std::cout << "Mutation rate: " << mutationRate << std::endl;
                    std::cout << "Random number: " << r << std::endl;
                }
                if (r < mutationRate) {
                    index = rand() % numItems;
                    if (debug){
                        std::cout << "Mutating index " << index << " of offspring ";
                        offspring[i].print();
                    }
                    offspring[i].chromosome.flip(index);
                }
                if (weight(offspring[i].chromosome) > knapsackMax) {
                    std::vector<int> items;
                    for (int j = 0; j < numItems; j++) {
                        if (offspring[i].chromosome[j]) {
                            items.push_back(j);
                        }
                    }
                    while (weight(offspring[i].chromosome) > knapsackMax) {
                        index = rand() % items.size();
                        offspring[i].chromosome.set(items[index], 0);
                        items.erase(items.begin() + index);
                    }
                }
                offspring[i].value = fitness(offspring[i].chromosome);
                offspring[i].weight = weight(offspring[i].chromosome);
            }
        };
        void replace(){
            // if (debug)
            //     std::cout << "Replacing" << std::endl;
            for (int i = 0; i < offspring.size(); i++) {
                bool chromosomeExists = exists(offspring[i].chromosome, &population);
                if (chromosomeExists == false && offspring[i].weight <= knapsackMax) {
                    population.push_back(offspring[i]);
                }
            }
            // if (debug)
            //     std::cout << "Population size before: " << population.size() << std::endl;
            std::sort(population.begin(), population.end(), 
                [](const Solution<111>& a, const Solution<111>& b) {
                    return a.value > b.value;
                }
            );
            double annoyance = 0.0;
            while (population.size() > populationSize) {
                //select random index above elitismCount
                // if (debug)
                //     std::cout << "Population size: " << population.size() << std::endl;
                int index = rand() % (population.size() - elitismCount) + elitismCount;
                // if (debug){
                //     std::cout << "Removing individual " << index << ": ";
                //     population[index].print();
                // }
                // remove with fitness proportionate probability
                double r = (double)rand() / RAND_MAX;
                if (r > population[index].value / (population[0].value * (1 + annoyance))) {
                    // if (debug){
                    //     std::cout << "Stochastic threshold met" << std::endl;
                    // }
                    population.erase(population.begin() + index);
                } 
                else {
                    // if (debug){
                    //     std::cout << r << " <= " << population[index].value / population[0].value << ". ";
                    //     std::cout << "Not removing: ";
                    //     population[index].print();
                    // }
                    annoyance += 0.05;
                }
            }
        };
        // greedy hill-climbing local search
        class LSSol {
            public:
                std::vector<int> items;
                double value = 0;
                double weight = 0;
                LSSol& operator=(const LSSol& sol) {
                    if (this == &sol) {
                        return *this;
                    }
                    this->items = sol.items;
                    this->value = sol.value;
                    this->weight = sol.weight;
                    return *this;
                }
        };
        void localSearch(Solution<111> &sol) {
            if (debug){
                std::cout << "Local search" << std::endl;
                std::cout << "Before: ";
                sol.print();
            }
            std::bitset<111> best = sol.chromosome;
            std::vector<int> solItems;
            std::stack<LSSol> searchSpace;
            // for (int i = 0; i < solItems.size(); i++) {
            //     searchSpace.push(solItems[i]);
            // }
            for (int i = 0; i < sol.length; i++) {
                if (sol.chromosome[i]) {
                    solItems.push_back(i);
                }
            }
            // sort by efficiency (value/weight) so we can pop the least efficient items
            LSSol currSol;
            LSSol bestSol;
            std::vector<LSSol> children;
            bestSol.items = solItems;
            bestSol.value = sol.value;
            bestSol.weight = sol.weight;
            searchSpace.push(bestSol);

            int i = 0;
            bool betterChildFound;
            do {
                betterChildFound = false;
                currSol = bestSol;
                // std::cout << "Current Solution: ";
                // for (int i = 0; i < currSol.items.size(); i++) {
                //     std::cout << currSol.items[i] << " ";
                // }
                // std::cout << " [V: " << currSol.value << "]";
                // std::cout << "[W: " << currSol.weight << "]\n";
                // searchSpace.pop();
                // children.clear();
                for (int i = 0; i < currSol.items.size(); ++i) {
                    LSSol child;
                    child = currSol;
                    // std::cout << "\t Items set\n";
                    child.items.erase(child.items.begin() + i);
                    // std::cout << "Item removed\n";
                    child.value = fitness(child.items);
                    child.weight = weight(child.items);
                    // add most efficient items
                    // std::cout << "\tChild " << i << ": ";
                    // for (int i = 0; i < child.items.size(); i++) {
                    //     std::cout << child.items[i] << " ";
                    // }
                    // std::cout << " [V: " << child.value << "]";
                    // std::cout << "[W: " << child.weight << "]\n";
                    for (int x = 0; x < numItems; ++x) {
                        if (std::find(child.items.begin(), child.items.end(), itemsEfficiencyIndex[x]) == child.items.end()
                            && child.weight + items[itemsEfficiencyIndex[x]][1] <= knapsackMax) 
                        {
                            child.items.push_back(itemsEfficiencyIndex[x]);
                            child.value = fitness(child.items);
                            child.weight = weight(child.items);
                        }
                    }
                    // std::cout << "\tChild " << i << " after: ";
                    // std::sort(child.items.begin(), child.items.end());
                    // for (int i = 0; i < child.items.size(); i++) {
                    //     std::cout << child.items[i] << " ";
                    // }
                    // std::cout << " [V: " << child.value << "]";
                    // std::cout << "[W: " << child.weight << "]\n";
                    if (child.value > bestSol.value) {
                        bestSol = child;
                        betterChildFound = true;
                    }
                }
                // if (betterChildFound) {
                //     std::cout << "Child " << i << " is better than parent" << std::endl;
                // }
            } while(betterChildFound);
            if (bestSol.value > sol.value) {
                sol.setChromosome(std::bitset<111>());
                for (int i = 0; i < bestSol.items.size(); i++) {
                    sol.chromosome.set(bestSol.items[i], 1);
                }
                sol.value = bestSol.value;
                sol.weight = bestSol.weight;
                if (debug) {
                    std::cout << "Improved solution: ";
                    sol.print();
                }
            }
        }
        void print(){
            for (int i = 0; i < populationSize; i++) {
                std::cout << i << ": ";
                // for (int j = 0; j < population[i].length; j++) {
                //     std::cout << "[" << items[population[i].chromosome[j]][0] << "|" << items[population[i].chromosome[j]][1] << "] ";
                // }
                // std::cout << ", Value: " << population[i].value;
                // std::cout << ", Weight: " << population[i].weight << std::endl;
                population[i].print();
            }
        };
};


std::unordered_map<
    std::string, 
    std::pair<
        std::vector<std::array<double, 3>>, 
        double
    >
> datasets;
class SearchConfig {
    public:
        int populationSize;
        int generations;
        double crossoverRate;
        double mutationRate;
        int elitismCount;
        std::string filename;
        double optimum;
        SearchConfig(
            int populationSize, 
            int generations, 
            double crossoverRate, 
            double mutationRate, 
            int elitismCount, 
            std::string filename,
            double optimum
        ) : populationSize(populationSize), 
            generations(generations), 
            crossoverRate(crossoverRate), 
            mutationRate(mutationRate), 
            elitismCount(elitismCount), 
            filename(filename),
            optimum(optimum)
        {}
};
int main() {
    DIR *dir;
    struct dirent *ent;
    std::string dirName = "./instances";
    if ((dir = opendir (dirName.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string fileName = ent->d_name;
            // Skip the file "graphs.xlsx" and any hidden files/directories
            if (fileName == "Known Optimums.xlsx" || fileName[0] == '.') {
                continue;
            }
            
            std::string filePath = dirName + "/" + fileName;
            std::ifstream file(filePath);
            if (file.is_open()) {
                std::cout << "Reading file: " << fileName << std::endl;
                std::string line;
                // get first line
                getline(file, line);
                std::istringstream iss(line);
                // first elem of line is number of remaining lines
                // second elem of line is max value
                int numLines;
                double max;

                if (!(iss >> numLines >> max)) {
                    std::cout << "Error reading first line of file: " << fileName << std::endl;
                    return EXIT_FAILURE;
                }

                datasets[fileName] = std::make_pair(std::vector<std::array<double, 3>>(), max);

                while (getline(file, line)) {
                    std::istringstream iss(line);
                    double x, y;
                    if (!(iss >> x >> y)) {
                        std::cout << "Error reading line of file: " << fileName << std::endl;
                        continue;
                    }
                    std::array<double, 3> point = {x, y, x/y};
                    datasets[fileName].first.push_back(point);
                }
                file.close();
                // // sort vector by y value
                std::sort(datasets[fileName].first.begin(), datasets[fileName].first.end(), 
                    [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
                        return a[1] < b[1];
                    }
                );
            } else {
                std::cout << "Unable to open file: " << fileName << std::endl;
            }
        }
        closedir (dir);
    } else {
        perror ("Unable to open directory");
        return EXIT_FAILURE;
    }

    // iterate over datasets
    for (auto dataset : datasets) {
        std::cout << "Dataset: " << dataset.first << std::endl;
        std::cout << "Max: " << dataset.second.second << ", num elem: " << dataset.second.first.size() << std::endl;
        for (std::array<double, 3> point : dataset.second.first) {
            std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
        }
    }

    std::vector<SearchConfig> searchConfigs({
        SearchConfig(10, 500, 0.9, 0.05, 1, "f1_l-d_kp_10_269", 295),
        SearchConfig(20, 500, 0.9, 0.05, 2, "f2_l-d_kp_20_878", 1024),
        SearchConfig(4, 500, 0.9, 0.05, 1, "f3_l-d_kp_4_20", 35),
        SearchConfig(4, 500, 0.9, 0.05, 1, "f4_l-d_kp_4_11", 23),
        SearchConfig(15, 500, 0.9, 0.05, 2, "f5_l-d_kp_15_375", 481.0694),
        SearchConfig(10, 500, 0.9, 0.05, 1, "f6_l-d_kp_10_60", 52),
        SearchConfig(7, 500, 0.9, 0.05, 1, "f7_l-d_kp_7_50",107),
        SearchConfig(23, 500, 0.9, 0.05, 3, "f8_l-d_kp_23_10000", 9767),
        SearchConfig(5, 500, 0.9, 0.05, 1, "f9_l-d_kp_5_80", 130),
        SearchConfig(20, 500, 0.9, 0.05, 2, "f10_l-d_kp_20_879", 1025),
        SearchConfig(100, 500, 0.9, 0.05, 10, "knapPI_1_100_1000_1", 9147)
    });
    bool runSeeded = false;
    if (runSeeded) {
        std::ofstream results("./results/seeded.csv");
        results << "Problem Instance,"
                << "Algorithm,"
                << "Seed Value,"
                << "Best Solution,"
                << "Known Optimum,"
                << "Runtime";
        for (SearchConfig config : searchConfigs) {
            double *avgBestPerGeneration = new double[config.generations];
            for (int i = 0; i < config.generations; i++) {
                avgBestPerGeneration[i] = 0;
            }
            double *avgBestPerGenerationMemetic = new double[config.generations];
            for (int i = 0; i < config.generations; i++) {
                avgBestPerGenerationMemetic[i] = 0;
            }
            double bestValue = 0;
            int seed = 0;

            for (bool memetic : {false, true}) {
                std::cout << "Running " << config.filename << (memetic ? " with memetic" : " without memetic") << std::endl;
                // generations << "Generation,Average Best\n";
                clock_t start = clock();
                GA ga(
                    datasets[config.filename].first, // items
                    datasets[config.filename].first.size(), // length
                    datasets[config.filename].second, // knapsackMax
                    config.populationSize, // populationSize
                    config.generations, // generations
                    config.crossoverRate, // crossoverRate
                    config.mutationRate, // mutationRate
                    config.elitismCount, // elitismCount
                    memetic, // memetic
                    false,
                    config.filename // filename (for seeding random number generator)
                );
                std::cout << "\tAverage weight: " << ga.avgWeight << std::endl;
                std::cout << "\tNum Items: " << ga.numItems << std::endl;
                std::cout << "\tPopulation size: " << ga.populationSize << std::endl;
                int solSize = round(ga.knapsackMax / ga.avgWeight);
                std::cout << "\tAverage solution size: " <<  solSize << "(" << round(ga.knapsackMax / ga.avgWeight) << ")" << std::endl;
                // roof of knapsackMax / avgWeight
                long long combinations = ga.numItems;
                for (int x = 1; x < solSize; x++) {
                    // std::cout << "\t\tCombination " << x << ": " << combinations << std::endl;
                    combinations *= ga.numItems - x;
                    combinations /= x + 1;
                }
                std::cout << "\tPotential number of solutions: " << combinations << std::endl;
                ga.run();
                double ratio = combinations / ga.populationSize;
                // std::cout << "\tMutation Rate: " << (0.01 + (0.1 - 0.01) * (ratio - 1) / (combinations - 1)) << std::endl;
                double logRatio = log10(ratio);
                double mutationRate = 0.01 + (0.1 - 0.01) * (logRatio - 0) / (log10(combinations) - 0);
                std::cout << "\tMutation Rate: " << mutationRate << std::endl;
                std::cout << "\tBest Sol: ";
                ga.population[0].print();
                bestValue = ga.population[0].value;
                clock_t end = clock();
                double time = (double)(end - start) / CLOCKS_PER_SEC;
                results << config.filename << ","
                        << (memetic ? "Memetic" : "Genetic") << ","
                        << ga.seed << ","
                        << bestValue << ","
                        << config.optimum << ","
                        << time << std::endl;
            }
        }
        results.close();
    }
    bool runAll = true;
    if (runAll)
    for (SearchConfig config : searchConfigs) {
        double *avgBestPerGeneration = new double[config.generations];
        for (int i = 0; i < config.generations; i++) {
            avgBestPerGeneration[i] = 0;
        }
        double *avgBestPerGenerationMemetic = new double[config.generations];
        for (int i = 0; i < config.generations; i++) {
            avgBestPerGenerationMemetic[i] = 0;
        }
        double bestValueNormal = 0;
        double avgBestNormal = 0;
        double bestValueMemetic = 0;
        double avgBestMemetic = 0;
        double avgTimeNormal = 0;
        double avgTimeMemetic = 0;
        double avgBestGenNormal = 0;
        double avgBestGenMemetic = 0;
        std::ofstream generations("./results/" + config.filename + "_generations.csv");
        generations << "Generation,Normal Average,Memetic Average\n";
        std::ofstream results("./results/" + config.filename + "_results.csv");
        results << "Average Time Normal," 
                << "Average Time Memetic," 
                << "Best Value Normal," 
                << "Best Value Memetic," 
                << "Average Best Normal," 
                << "Average Best Memetic," 
                << "Generations," 
                << "Population Size," 
                << "Crossover Rate," 
                << "Mutation Rate," 
                << "Elitism Count," 
                << "Average Generation of Best Value Normal," 
                << "Average Generation of Best Value Memetic\n";

        int numIterations = 40;
        for (bool memetic : {false, true}) {
            std::cout << "Running " << config.filename << (memetic ? " with memetic" : " without memetic") << std::endl;
            // generations << "Generation,Average Best\n";
            for (int i = 0; i < numIterations; ++i) {
                clock_t start = clock();
                GA ga(
                    datasets[config.filename].first, // items
                    datasets[config.filename].first.size(), // length
                    datasets[config.filename].second, // knapsackMax
                    config.populationSize, // populationSize
                    config.generations, // generations
                    config.crossoverRate, // crossoverRate
                    config.mutationRate, // mutationRate
                    config.elitismCount, // elitismCount
                    memetic, // memetic
                    false
                    // config.filename // filename (for seeding random number generator)
                );
                ga.run();
                clock_t end = clock();
                config.mutationRate = ga.mutationRate;
                config.crossoverRate = ga.crossoverRate;
                config.populationSize = ga.populationSize;
                config.generations = ga.generations;
                double time = (double)(end - start) / CLOCKS_PER_SEC;
                if (memetic) {
                    avgBestMemetic += ga.population[0].value;
                    avgBestGenMemetic += ga.generationOfBest;
                    avgTimeMemetic += time;
                    if (ga.population[0].value > bestValueMemetic) {
                        bestValueMemetic = ga.population[0].value;
                    }
                } else {
                    avgBestNormal += ga.population[0].value;
                    avgBestGenNormal += ga.generationOfBest;
                    avgTimeNormal += time;
                    if (ga.population[0].value > bestValueNormal) {
                        bestValueNormal = ga.population[0].value;
                    }
                }
                for (int i = 0; i < config.generations; i++) {
                    if (memetic) {
                        avgBestPerGenerationMemetic[i] += ga.bestFitnessPerGeneration[i];
                    } else {
                        avgBestPerGeneration[i] += ga.bestFitnessPerGeneration[i];
                    }
                }
            }
        }
        avgTimeNormal /= numIterations;
        avgTimeMemetic /= numIterations;
        avgBestNormal /= numIterations;
        avgBestMemetic /= numIterations;
        avgBestGenNormal /= numIterations;
        avgBestGenMemetic /= numIterations;
        results << avgTimeNormal << ","
                << avgTimeMemetic << "," 
                << bestValueNormal << "," 
                << bestValueMemetic << "," 
                << avgBestNormal << "," 
                << avgBestMemetic << "," 
                << config.generations << "," 
                << config.populationSize << "," 
                << config.crossoverRate << "," 
                << config.mutationRate << "," 
                << config.elitismCount << "," 
                << avgBestGenNormal << "," 
                << avgBestGenMemetic << "\n";
        for (int i = 0; i < config.generations; i++) {
            double normal = avgBestPerGeneration[i] / (double)numIterations;
            double memetic = avgBestPerGenerationMemetic[i] / (double)numIterations;
            // std::cout << "\t" << i << ": " << memetic << std::endl;
            generations << i << "," << normal << "," << memetic << "\n";
        }
        delete [] avgBestPerGeneration;
        delete [] avgBestPerGenerationMemetic;
        results.close();
        generations.close();
    }
    return 0;
}
