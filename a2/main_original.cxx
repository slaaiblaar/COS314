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
class Solution {
    public:
        Solution(const Solution& sol): i(sol.i), length(sol.length), value(sol.value), weight(sol.weight) {
            if (sol.chromosome != nullptr) {
                this->chromosome = new int[sol.length];
                for (int i = 0; i < sol.length; i++) {
                    this->chromosome[i] = sol.chromosome[i];
                }
            }
        }
        Solution& operator=(const Solution& sol) {
            if (this == &sol) {
                return *this;
            }
            this->i = sol.i;
            this->length = sol.length;
            this->value = sol.value;
            this->weight = sol.weight;
            if (this->chromosome != nullptr) {
                delete[] this->chromosome;
                this->chromosome = nullptr;
            }
            if (sol.chromosome != nullptr) {
                this->chromosome = new int[sol.length];
                for (int i = 0; i < sol.length; i++) {
                    this->chromosome[i] = sol.chromosome[i];
                }
            }
            // sort
            std::sort(this->chromosome, this->chromosome + this->length);
            return *this;
        }
        Solution(){};
        ~Solution() {
            if (this->chromosome != nullptr)
            {
                delete[] this->chromosome;
                this->chromosome = nullptr;
            }
        }
        void setChromosome(std::vector<int> chromosome) {
            std::sort(chromosome.begin(), chromosome.end());
            if (this->chromosome != nullptr)
            {
                delete [] this->chromosome;
                this->chromosome = nullptr;
            }
            this->chromosome = new int[chromosome.size()];
            for (int i = 0; i < chromosome.size(); i++) {
                this->chromosome[i] = chromosome[i];
            }
            this->length = chromosome.size();
            this->value = 0;
            this->weight = 0;
        }
        void setChromosome(int *chromosome, int length) {
            std::sort(chromosome, chromosome + length);
            if (this->chromosome != nullptr)
            {
                delete [] this->chromosome;
                this->chromosome = nullptr;
            }
            this->chromosome = new int[length];
            for (int i = 0; i < length; i++) {
                this->chromosome[i] = chromosome[i];
            }
            this->length = length;
            this->value = 0;
            this->weight = 0;
        }
        void print() {
            for (int i = 0; i < length; i++) {
                std::cout << chromosome[i] << " ";
            }
            std::cout << " [V: " << value << "]";
            std::cout << "[W: " << weight << "]\n";
        }
        int i;
        int *chromosome = nullptr;
        int length = 0;
        double value;
        double weight;
};
class GA {
    public:
        std::string filename;
        bool debug = false;
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
            std::string filename = ""
        ) : filename(filename),
            numItems(length), 
            knapsackMax(knapsackMax), 
            populationSize(populationSize), 
            generations(generations), 
            crossoverRate(crossoverRate), 
            mutationRate(mutationRate), 
            elitismCount(elitismCount),
            memetic(memetic)
        {
            // debug = (filename == "f6_l-d_kp_10_60" && populationSize == 15 && generations == 100 && crossoverRate == 0.5 && mutationRate == 0.07 && elitismCount == 7);
            debug = false;
            unsigned int seed = 0;
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
            }
            // std::cout << "Items populated" << std::endl;
            this->init();
            // std::cout << "Init complete" << std::endl;
        }
        ~GA(){
            // std::cout << "GA destructor" << std::endl;
            delete [] items;
            // delete [] population;

        }
        void run(){
            double bestFitness = 0;
            for (int i = 0; i < generations; i++) {
                // std::cout << "Generation " << i << std::endl;
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
                replace();
                if (debug)
                    std::cout << "Replace complete" << std::endl;
            }
            // print();
            // std::cout << "Final population: " << std::endl;
            // for (int i = 0; i < populationSize; i++) {
            //     std::cout << "\t" << i << ": ";
            //     population[i].print();
            // }
            std::cout << "Gen " << std::setw(4) << std::setfill('0') << generationOfBest << ". Best sol: ";
            population[0].print();
            
            // std::cout << "Done" << std::endl;
        }
        std::vector<Solution> population;
        std::vector<Solution> offspring;
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
        double totalFitness;
        int generationOfBest;

    private:
        // pair->first is chromosome, pair->second is fitness
        // Solution *population;
        void init() {
            // std::cout << "Initialising population" << std::endl;
            // std::cout << "numItems: " << numItems << std::endl;
            // std::cout << "knapsackMax: " << knapsackMax << std::endl;
            for (int i = 0; i < populationSize; i++) {
                bool *present = new bool[numItems];
                for (int j = 0; j < numItems; j++) {
                    present[j] = false;
                }
                std::vector<int> chromosomeVec;
                // std::cout << "Chromosome " << i << std::endl;
                // std::cout << "numItems: " << numItems << std::endl;
                // std::cout << "knapsackMax: " << knapsackMax << std::endl;
                // random chromosome length grater than 2
                int chromosomeLength = rand() % (numItems - 3) + 3;
                // std::cout << "Generating chromosome of length " << chromosomeLength << std::endl;
                for (int j = 0; j < numItems
                 && chromosomeVec.size() < chromosomeLength
                 && weight(chromosomeVec) < knapsackMax; j++) {
                    // std::cout << "Here 1\n";
                    int index = rand() % numItems;
                    int counter = 0;
                    while (index > 0 && counter < numItems) {
                        if (!present[counter]) {
                            --index;
                        }
                        ++counter;
                    }
                    // std::cout << "Counter: " << counter << std::endl;
                    while (counter < numItems && present[counter]) {
                        counter++;
                    }
                    if (counter == numItems) {
                        continue;
                    }
                    // std::cout << "Here 2\n";
                    present[counter] = true;
                    // std::cout << "Here 3\n";
                    chromosomeVec.push_back(counter);
                }
                // std::cout << "Here 4\n";
                while (chromosomeVec.size() > 0 && weight(chromosomeVec) > knapsackMax) {
                    chromosomeVec.pop_back();
                }
                // std::cout << "Here 5\n";
                if (chromosomeVec.size() == 0) {
                    // std::cout << "Chromosome " << i << " has no items" << std::endl;
                    i--;
                    delete[] present;
                    continue;
                }
                // std::cout << "Here 6\n";
                double weight = this->weight(chromosomeVec);

                if (weight == 0) {
                    // std::cout << "Chromosome " << i << " has weight 0" << std::endl;
                    i--;
                    delete[] present;
                    continue;
                }
                // std::cout << "Here 7\n";
                bool chromosomeExists = exists(chromosomeVec);
                if (chromosomeExists) {
                    // std::cout << "Chromosome " << i << " already exists" << std::endl;
                    i--;
                    delete[] present;
                    continue;
                }
                Solution sol;
                sol.setChromosome(chromosomeVec);
                sol.weight = weight;
                sol.value = fitness(chromosomeVec);
                population.push_back(sol);
                delete[] present;
            }
            
            std::sort(population.begin(), population.end(), 
                [](const Solution& a, const Solution& b) {
                    return a.value > b.value;
                }
            );
            // std::cout << "Initial population: \n";
            // for (int i = 0; i < population.size(); ++i) {
            //     population[i].print();
            // }
        }
        bool exists(std::vector<int> chromosome, std::vector<Solution> *pop = nullptr) {
            if (pop == nullptr) {
                pop = &population;
            }
            std::sort(chromosome.begin(), chromosome.end());
            for (int i = 0; i < pop->size() && (*pop)[i].length != 0; i++) {
                bool found = true;
                if ((*pop)[i].length != chromosome.size()
                //  || fitness((*pop)[i].chromosome, (*pop)[i].length) != fitness(chromosome)
                //  || weight((*pop)[i].chromosome, (*pop)[i].length) != weight(chromosome)
                ) {
                    continue;
                }
                for (int j = 0; j < chromosome.size(); j++) {
                    if ((*pop)[i].chromosome[j] != chromosome[j]) {
                        break;
                    }
                }
                return true;
            }
            return false;
        }
        double fitness(int *chromosome, int length) {
            double value = 0;
            for (int i = 0; i < length; i++) {
                value += items[chromosome[i]][0];
            }
            return value;
        }
        double fitness(std::vector<int> chromosome) {
            double value = 0;
            for (int i = 0; i < chromosome.size(); i++) {
                value += items[chromosome[i]][0];
            }
            return value;
        }
        double weight(int *chromosome, int length) {
            double weight = 0;
            for (int i = 0; i < length; i++) {
                weight += items[chromosome[i]][1];
            }
            return weight;
        }
        double weight(std::vector<int> chromosome) {
            double weight = 0;
            for (int i = 0; i < chromosome.size(); i++) {
                weight += items[chromosome[i]][1];
            }
            return weight;
        }

        void eval(){
            totalFitness = 0;
            for (int i = 0; i < populationSize; i++) {
                totalFitness += population[i].value;
            }
            std::sort(population.begin(), population.end(), 
                [](const Solution& a, const Solution& b) {
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
            int numSelected = 0;
            while (numSelected < numParents) {
                int index = rand() % populationSize;
                if (!selected[index]) {
                    double r = (double)rand() / RAND_MAX;
                    if (r < population[index].value / totalFitness) {
                        parents.push_back(index);
                        selected[index] = true;
                        numSelected++;
                    }
                }
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
                    // std::cout << "Crossover threshold met" << std::endl;
                    int crossoverPoint = rand() % population[currParents[0]].length;
                    // std::cout << "Crossover point: " << crossoverPoint << std::endl;
                    std::vector<int> child1;
                    std::vector<int> child2;
                    for (int i = 0; i < crossoverPoint; i++) {
                        // add if not present
                        if (std::find(child1.begin(), child1.end(), population[currParents[0]].chromosome[i]) == child1.end()) {
                            child1.push_back(population[currParents[0]].chromosome[i]);
                        }
                        if (std::find(child2.begin(), child2.end(), population[currParents[1]].chromosome[i]) == child2.end()) {
                            child2.push_back(population[currParents[1]].chromosome[i]);
                        }
                    }
                    // std::cout << "Child 1: ";
                    // for (int i = 0; i < child1.size(); i++) {
                    //     std::cout << child1[i] << " ";
                    // }
                    // std::cout << std::endl;
                    // std::cout << "Child 2: ";
                    // for (int i = 0; i < child2.size(); i++) {
                    //     std::cout << child2[i] << " ";
                    // }
                    // std::cout << std::endl;
                    for (int i = crossoverPoint; i < population[currParents[0]].length; i++) {
                        // add if not present
                        if (std::find(child2.begin(), child2.end(), population[currParents[0]].chromosome[i]) == child2.end()) {
                            child2.push_back(population[currParents[0]].chromosome[i]);
                        }
                        if (std::find(child1.begin(), child1.end(), population[currParents[1]].chromosome[i]) == child1.end()) {
                            child1.push_back(population[currParents[1]].chromosome[i]);
                        }
                    }
                    while (fitness(child1) > knapsackMax) {
                        // remove random item
                        int index = rand() % child1.size();
                        child1.erase(child1.begin() + index);
                    }
                    while (fitness(child2) > knapsackMax) {
                        // remove random item
                        int index = rand() % child2.size();
                        child2.erase(child2.begin() + index);
                    }
                    if (child1.size() > 0) {
                        Solution sol1;
                        sol1.setChromosome(child1);
                        sol1.value = fitness(child1);
                        sol1.weight = weight(child1);
                        offspring.push_back(sol1);
                    }
                    if (child2.size() > 0) {
                        Solution sol2;
                        sol2.setChromosome(child2);
                        sol2.value = fitness(child2);
                        sol2.weight = weight(child2);
                        offspring.push_back(sol2);
                    }
                } 
                // else 
                // {
                //     // as per Sastry, K., Goldberg, D., & Kendall, G. (2014). Genetic algorithms. In Search methodologies (pp. 96). Springer, Boston, MA.
                //     Solution sol1;
                //     Solution sol2;
                //     sol1 = population[currParents[0]];
                //     sol2 = population[currParents[1]];
                //     offspring.push_back(sol1);
                //     offspring.push_back(sol2);
                // }
            }
            // std::cout << "Offspring: " << std::endl;
            // for (int i = 0; i < offspring.size(); i++) {
            //     std::cout << "\tO " << i << ": ";
            //     for (int j = 0; j < offspring[i].length; j++) {
            //         std::cout << offspring[i].chromosome[j] << " ";
            //     }
            //     std::cout << ", Value: " << offspring[i].value;
            //     std::cout << ", Weight: " << offspring[i].weight << std::endl;
            // }
        }
        void mutate(){
            // Possible Mutations: Adding allele, Removing Allele, Each with separate probabilities
            bool selected[numItems];
            int numTested = 0;
            int index;
            double r;
            for (int i = 0; i < offspring.size(); i++) {
                if (debug){
                    std::cout << "Offspring " << i << std::endl;
                    offspring[i].print();
                }
                numTested = 0;
                r = (double)rand() / RAND_MAX;
                if (r < mutationRate) {
                    // std::cout << "Mutating offspring: ";
                    // offspring[i].print();
                    std::vector<int> temp(offspring[i].chromosome, offspring[i].chromosome + offspring[i].length);
                    // mutate
                    for (int j = 0; j < numItems; j++) {
                        selected[j] = false;
                    }
                    for (int j = 0; j < offspring[i].length; j++) {
                        selected[offspring[i].chromosome[j]] = true;
                        ++numTested;
                    }
                    // for (int j = 0; j < numItems; j++) {
                    //     std::cout << "Item " << j << " selected: " << selected[j] << std::endl;
                    // }
                    // std::cout << "Num items tested: " << numTested << std::endl;
                    // std::cout << "Num items: " << numItems << std::endl;
                    int mutationType = rand() % 3;
                    // std::cout << "Mutation type: " << mutationType << std::endl;
                    switch (mutationType) {
                        case 0:
                            // add allele
                            do {
                                int counter;
                                do {
                                    counter = 0;
                                    int numAvailable = numItems - numTested;
                                    // std::cout << "Num available: " << numAvailable << std::endl;
                                    index = rand() % numAvailable;
                                    // std::cout << "Index before: " << index << std::endl;
                                    while (index > 0 && counter < numItems) {
                                        if (!selected[counter]) {
                                            --index;
                                        }
                                        ++counter;
                                    }
                                    // std::cout << counter << std::endl;
                                    while (counter < numItems && selected[counter]) {
                                        ++counter;
                                    }
                                    // std::cout << "Counter: " << counter << std::endl;
                                    if (counter == numItems) {
                                        continue;
                                    }
                                    // std::cout << "index after: " << index << std::endl;
                                    // std::cout << "Counter " << counter << " selected: " << selected[counter] << std::endl;

                                } while (selected[counter] && numTested < numItems);
                                selected[index] = true;
                                ++numTested;
                                index = counter;
                            } while (weight(temp) + items[index][1] > knapsackMax && numTested < numItems);
                            if (numTested == numItems) {
                                break;
                            }
                            temp.push_back(index);
                            break;
                        case 1:
                            // remove allele
                            index = rand() % offspring[i].length;
                            temp.erase(temp.begin() + index);
                            break;
                        default:
                            break;
                    }
                    while (weight(temp) > knapsackMax) {
                        // remove random item
                        int index = rand() % temp.size();
                        temp.erase(temp.begin() + index);
                    }
                    offspring[i].setChromosome(temp);
                    offspring[i].value = fitness(temp);
                    offspring[i].weight = weight(temp);
                    // std::cout << "Mutated offspring: ";
                    // offspring[i].print();
                }
            }
        };
        void replace(){
            // if (debug)
            //     std::cout << "Replacing" << std::endl;
            for (int i = 0; i < offspring.size(); i++) {
                std::vector<int> chromosomeVec(offspring[i].chromosome, offspring[i].chromosome + offspring[i].length);
                bool chromosomeExists = exists(chromosomeVec);
                if (chromosomeExists == false) {
                    population.push_back(offspring[i]);
                }
            }
            // if (debug)
            //     std::cout << "Population size before: " << population.size() << std::endl;
            std::sort(population.begin(), population.end(), 
                [](const Solution& a, const Solution& b) {
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
        void selectNextGen(){};
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


    double avgTime = 0;
    for (int i = 0; i < 20; ++i) {
        clock_t start = clock();
        // GA(items, length, knapsackMax, populationSize, generations, crossoverRate, mutationRate, elitismCount, memetic, filename)
        // GA ga(
        //     datasets["f5_l-d_kp_15_375"].first, datasets["f5_l-d_kp_15_375"].first.size(), datasets["f5_l-d_kp_15_375"].second, 
        //     45, 1000, 0.9, 0.1, 4, false
        //     // "f5_l-d_kp_15_375" // filename (for seeding random number generator)
        // );
        GA ga(
            datasets["knapPI_1_100_1000_1"].first, datasets["knapPI_1_100_1000_1"].first.size(), datasets["knapPI_1_100_1000_1"].second, 
            100, 100, 0.5, 0.05, 10, false
            // "f5_l-d_kp_15_375" // filename (for seeding random number generator)
        );
        ga.run();
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        avgTime += time;
    }
    avgTime /= 20;
    std::cout << "Average time taken: " << avgTime << "s" << std::endl;

    // Run GA on each dataset with many different configurations as follows:
    // 1. Population Sizes between 5 and thrice the number of items in the dataset
    // 2. Generations between 100 and 1000
    // 3. Crossover Rates between 0.5 and 0.9
    // 4. Mutation Rates between 0.01 and 0.1
    // 5. Elitism Counts between 1 and half the population size
    // 6. Memetic Algorithm on and off
    // Store the following in a csv file named "${filename}_results.csv" and "${filename}_results_memetic.csv":
    // 1. Average time
    // 2. Best value
    // 3. Number of generations
    // 4. Population Size
    // 5. Crossover Rate
    // 6. Mutation Rate
    // 7. Elitism Count
    // 8. Generation in which best value was found
    if (true)
    for (auto dataset : datasets) {
        std::string filename = dataset.first;
        std::ofstream results
        ("./results/" + filename + "_results.csv");
        results << "Average Time,Best Value,Generations,Population Size,Crossover Rate,Mutation Rate,Elitism Count,Generation of Best Value\n";
        std::ofstream resultsMemetic
        ("./results/" + filename + "_results_memetic.csv");
        resultsMemetic << "Average Time,Best Value,Generations,Population Size,Crossover Rate,Mutation Rate,Elitism Count,Generation of Best Value\n";
        for (int populationSize = dataset.second.first.size(); populationSize <= dataset.second.first.size() * 3; populationSize += dataset.second.first.size()) {
            for (int generations = 100; generations <= 300; generations += 100) {
                for (double crossoverRate = 0.5; crossoverRate <= 0.9; crossoverRate += 0.1) {
                    for (double mutationRate = 0.05; mutationRate <= 0.1; mutationRate += 0.01) {
                        // for (int elitismCount = 1; elitismCount <= populationSize / 2; elitismCount++) {
                            // elitism is max(1, 0.1*populationSize)
                            int elitismCount = std::max(1, (int) (0.1 * populationSize));
                            // for (bool memetic : {true, false}) {
                                bool memetic = false;
                                double avgTime = 0;
                                double bestValue = 0;
                                int bestGeneration = 0;
                                std::cout << "Current configuration: file: " << filename << "pop size: " << populationSize << ", gens: " << generations << ", x rate: " << crossoverRate << ", m rate: " << mutationRate << ", elites: " << elitismCount << ", memetic: " << memetic << std::endl;
                                for (int i = 0; i < 20; ++i) {
                                    clock_t start = clock();
                                    GA ga(
                                        dataset.second.first, // items
                                        dataset.second.first.size(), // length
                                        dataset.second.second, // knapsackMax
                                        populationSize, // populationSize
                                        generations, // generations
                                        crossoverRate, // crossoverRate
                                        mutationRate, // mutationRate
                                        elitismCount, // elitismCount
                                        memetic // memetic
                                        // filename // filename (for seeding random number generator)
                                    );
                                    ga.run();
                                    clock_t end = clock();
                                    double time = (double)(end - start) / CLOCKS_PER_SEC;
                                    avgTime += time;
                                    if (ga.population[0].value > bestValue) {
                                        bestValue = ga.population[0].value;
                                        bestGeneration = ga.generationOfBest;
                                    }
                                }
                                avgTime /= 20;
                                if (memetic) {
                                    resultsMemetic << avgTime << "," << bestValue << "," << generations << "," << populationSize << "," << crossoverRate << "," << mutationRate << "," << elitismCount << "," << bestGeneration << "\n";
                                } else {
                                    results << avgTime << "," << bestValue << "," << generations << "," << populationSize << "," << crossoverRate << "," << mutationRate << "," << elitismCount << "," << bestGeneration << "\n";
                                }
                            // }
                        // }
                    }
                }
            }
        }
        results.close();
        resultsMemetic.close();
    }


    // clock_t start = clock();
    // GA ga(datasets["f5_l-d_kp_15_375"].first, datasets["f5_l-d_kp_15_375"].first.size(), datasets["f5_l-d_kp_15_375"].second, 500, 100, 0.8, 0.05, 10, false, "f5_l-d_kp_15_375");
    // // ga.init();
    // // for (int popI = 0; popI < ga.population.size(); popI++) {
    // //     std::cout << popI << ": ";
    // //     for (int i = 0; i < ga.population[popI].length; i++) {
    // //         std::cout << ga.population[popI].chromosome[i] << " ";
    // //     }
    // //     std::cout << ", Value: " << ga.population[popI].value;
    // //     std::cout << ", Weight: " << ga.population[popI].weight << std::endl;
    // // }
    // // start time
    // ga.run();
    // // end time
    // clock_t end = clock();
    // double time = (double)(end - start) / CLOCKS_PER_SEC;
    // std::cout << "Time taken: " << time << "s" << std::endl;

    return 0;
}
