#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <cstdlib>
#include <time.h>
#include <climits>
#include <cmath>
#include <sys/time.h>
#include <fstream>
#include <chrono>
// for swapping search spaces automatically
#ifdef MAX_ITERATIONS
const int MAX_CONSECUTIVE_ITERATIONS = MAX_ITERATIONS;
#else
const int MAX_CONSECUTIVE_ITERATIONS = 50;
#endif
#ifndef SEARCH_SPACE
#define SEARCH_SPACE CAMPUSES
#endif
#define CAMPUSES 5
#define LARGE_SPACE 29
#if SEARCH_SPACE == CAMPUSES
const int numVecs = 5;
int distances[5][5] = {
    {0, 15, 20, 22, 30},
    {15, 0, 10, 12, 25},
    {20, 10, 0, 8, 22},
    {22, 12, 8, 0, 18},
    {30, 25, 22, 18, 0}};
std::string labels[] = {"Hatfield", "Hillcrest", "Groenkloof", "Prinsof", "Mamelodi"};
#elif SEARCH_SPACE == LARGE_SPACE
const int numVecs = 29;
double distances[numVecs][numVecs] = {0};
double vectors[29][2] = {
    {519.5, 378},
    {371.5, 351},
    {210, 472.5},
    {405.5, 193},
    {405.5, 453.5},
    {485.5, 465.5},
    {656.5, 68.5},
    {614.5, 341.5},
    {417.5, 519},
    {393.5, 251.5},
    {434, 39.5},
    {526.5, 531},
    {468.5, 261},
    {337.5, 83},
    {405.5, 139},
    {558, 220},
    {259.5, 51.5},
    {321, 127},
    {490, 153.5},
    {359.5, 275.5},
    {429.5, 382.5},
    {330.5, 25},
    {706, 232},
    {551, 307},
    {558, 107.5},
    {330.5, 484.5},
    {605, 285},
    {551, 419.5},
    {294, 439}};
std::string labels[numVecs] = {""};
void calculateVecDists()
{
    for (int x = 0; x < 29; ++x)
    {
        labels[x] = std::to_string(x);
        labels[x] = std::string(2 - labels[x].length(), '0') + labels[x];
        for (int y = 0; y < 29; ++y)
        {
            distances[x][y] = sqrt(pow(vectors[x][0] - vectors[y][0], 2) + pow(vectors[x][1] - vectors[y][1], 2));
        }
    }
}
#endif
int heuristic(int *sol, int len)
{
    int totalVal = 0;
    for (int x = 0; x < numVecs - 1; ++x)
    {
        totalVal += distances[sol[x]][sol[x + 1]];
    }
    totalVal += distances[sol[numVecs - 1]][sol[0]];
    return totalVal;
}
// uses a greedy approach to swap any vertices that would produce a lower heuristic value
void optimiseLocally(int *sol, int len)
{
    int solHeur = heuristic(sol, len);
    int newHeur;
    int newSol[len];
    for (int x = 0; x < len; ++x)
    {
        newSol[x] = sol[x];
    }
    bool swapped = true;
    // Big O go BRRRR
    while (swapped)
    {
        swapped = false;
        for (int i = 1; i < len - 1; ++i)
        {
            newSol[i] = sol[i + 1];
            newSol[i + 1] = sol[i];
            newHeur = heuristic(newSol, len);
            // swap
            if (newHeur < solHeur)
            {
                sol[i] = newSol[i];
                sol[i + 1] = newSol[i + 1];
                solHeur = newHeur;
                if (i > 1)
                    --i;
            }
            else // reset
            {
                newHeur = solHeur;
                newSol[i] = sol[i];
                newSol[i + 1] = sol[i + 1];
            }
        }
    }
}
void perturb(int *sol, int len, int *perturbedSol)
{
    for (int x = 0; x < len; ++x)
    {
        perturbedSol[x] = sol[x];
    }
    int i = (rand() % (len - 1)) + 1;
    int j;
    do
    {
        j = (rand() % (len - 1)) + 1;
    } while (i - j < 1 && i - j > -1);

    perturbedSol[i] = sol[j];
    perturbedSol[j] = sol[i];
}

void generateInitialSol(int *sol)
{
    bool visited[numVecs] = {0};
    visited[0] = 1;
    int randCampus;
    for (int n = numVecs - 1; n > 0; --n)
    {
        randCampus = rand() % n;
        int counter = 0;
        while (randCampus >= 0 && counter < numVecs)
        {
            if (!visited[counter])
            {
                if (randCampus == 0)
                    break;
                --randCampus;
            }
            ++counter;
        }
        sol[numVecs - n] = counter;
        visited[counter] = true;
    }
}
int main()
{
    srand(time(NULL));

#if SEARCH_SPACE == LARGE_SPACE
    calculateVecDists();
#endif
    int currSol[numVecs] = {0};
    int padding[numVecs];
    int bestSol[numVecs];
    int padding2[numVecs];
    int consecutiveIterations;
    int padding3;
    timespec timeObj;
    std::chrono::_V2::system_clock::time_point time1; // disgusting
    std::chrono::_V2::system_clock::time_point time2; // disgusting
    int v1 = 0;
    int v2 = 0;
#ifdef RUN_MANY
    int iterations = 0;
    int solutions[RUN_MANY][numVecs + 2];
    int randomPadding[RUN_MANY][numVecs + 2];
    int avgPerturbedHeuristics[RUN_MANY][MAX_CONSECUTIVE_ITERATIONS];
    int moreRandomPadding[RUN_MANY][MAX_CONSECUTIVE_ITERATIONS];
    long int executionTime[RUN_MANY];
    long int evenMoreRandomPadding[RUN_MANY];
    int totalDiffCurrBest[RUN_MANY] = {0}; // Sum of (curHeur - bestHeur) for all iterations
    int diffPadding[RUN_MANY];
    int maxConsecutiveIterations = 1;
    double avgCostPerIterations[MAX_CONSECUTIVE_ITERATIONS] = {0};
    long double avgTimePerIterations[MAX_CONSECUTIVE_ITERATIONS] = {0};
    int bestCostPerIterations[MAX_CONSECUTIVE_ITERATIONS];
    for (int x = 0; x < MAX_CONSECUTIVE_ITERATIONS; ++x)
    {
        bestCostPerIterations[x] = INT_MAX;
    }

run_many:
    if (iterations >= RUN_MANY)
    {
        avgCostPerIterations[maxConsecutiveIterations - 1] /= (RUN_MANY * 1.0);
        avgTimePerIterations[maxConsecutiveIterations - 1] /= (RUN_MANY * 1.0);
        std::cout << maxConsecutiveIterations << ": " << avgCostPerIterations[maxConsecutiveIterations - 1] << " cost, " << (avgTimePerIterations[maxConsecutiveIterations - 1] / 1000.0) << "ms avg\n";
        std::cout << "Max iterations: " << maxConsecutiveIterations << std::endl;
        ++maxConsecutiveIterations;
        iterations = 0;
    }
    if (maxConsecutiveIterations > MAX_CONSECUTIVE_ITERATIONS)
        goto end_loop;
#else
    int maxConsecutiveIterations = MAX_CONSECUTIVE_ITERATIONS;
#endif
    time1 = std::chrono::high_resolution_clock::now();
    generateInitialSol(currSol);
    std::cout << "Initial Solution: ";
    for (int x = 0; x < numVecs; ++x)
    {
        std::cout << currSol[x] << " > ";
    }
    std::cout << currSol[0] << std::endl;
    std::cout << "Heuristic Value: " << heuristic(currSol, numVecs) << "\n";
    optimiseLocally(currSol, numVecs);
    std::cout << "Optimised initial solution: ";
    for (int x = 0; x < numVecs - 1; ++x)
    {
        bestSol[x] = currSol[x];
        std::cout << labels[bestSol[x]] << ">";
    }
    bestSol[numVecs - 1] = currSol[numVecs - 1];
    std::cout << labels[currSol[numVecs - 1]] << ">" << labels[0] << std::endl;
    std::cout << "Heuristic Value: " << heuristic(currSol, numVecs) << std::endl;
    consecutiveIterations = 0;
    while (consecutiveIterations <= maxConsecutiveIterations)
    {
        ++consecutiveIterations;
        for (int x = 1; x < numVecs; ++x)
        {
            currSol[x] = bestSol[x];
        }
        v1 = (rand() % (numVecs - 1)) + 1;
        do
        {
            v2 = (rand() % (numVecs - 1)) + 1;
        } while (v1 - v2 < 2 && v1 - v2 > -2);
        currSol[v1] = bestSol[v2];
        currSol[v2] = bestSol[v1];

        std::cout << "Perturbed Solution: \n\t";
        for (int x = 0; x < numVecs; ++x)
        {
            std::cout << currSol[x] << ">";
        }
        std::cout << currSol[0] << "\n\t";
        std::cout << "Heuristic: " << heuristic(currSol, numVecs) << std::endl;
        optimiseLocally(currSol, numVecs);

        std::cout << "Perturbed Solution After local Search: \n\t";
        for (int x = 0; x < numVecs; ++x)
        {
            std::cout << currSol[x] << ">";
        }
        std::cout << currSol[0] << "\n\t";
        std::cout << "Heuristic: " << heuristic(currSol, numVecs) << std::endl;
        if (heuristic(currSol, numVecs) < heuristic(bestSol, numVecs))
        {
            std::cout << "\033[42mNew optimal solution: \033[0m\n";
            for (int x = 0; x < numVecs; ++x)
            {
                std::cout << (bestSol[x] = currSol[x]) << ">";
            }
            std::cout << bestSol[0] << std::endl;
            std::cout << "Heuristic value: " << heuristic(bestSol, numVecs) << std::endl;
            consecutiveIterations = 0;
        }
    }
    time2 = std::chrono::high_resolution_clock::now();
#ifndef RUN_MANY
    std::cout << std::endl
              << "Optimal solution: \n";
    for (int x = 0; x < numVecs; ++x)
    {
        std::cout << labels[bestSol[x]] << ">";
    }

    std::cout << std::endl
              << labels[0]
              << "\nCost: " << heuristic(bestSol, numVecs) << std::endl;
#endif
#ifdef RUN_MANY // Ignore, just for data collection
    executionTime[iterations] = std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count();
    solutions[iterations][0] = heuristic(bestSol, numVecs);
    avgCostPerIterations[maxConsecutiveIterations - 1] += heuristic(bestSol, numVecs);
    avgTimePerIterations[maxConsecutiveIterations - 1] += executionTime[iterations];

    if (heuristic(bestSol, numVecs) < bestCostPerIterations[maxConsecutiveIterations - 1])
    {
        bestCostPerIterations[maxConsecutiveIterations - 1] = heuristic(bestSol, numVecs);
    }
    for (int x = 1; x <= numVecs; ++x)
    {
        solutions[iterations][x] = bestSol[x - 1];
    }
    solutions[iterations][numVecs + 1] = 0;
    ++iterations;
    goto run_many;
end_loop:
    double avgHeur = 0;
    double heurPadding = 0;
    long double avgTime = 0;
    long double timePadding = 0;
    int max[numVecs + 2] = {0};
    max[0] = INT16_MAX;
    std::ofstream csv;
    csv.open("ils_data.csv");
    csv << "Iteration, Time, Cost, Solution \n";
    for (int i = 0; i < RUN_MANY; ++i)
    {
        if (solutions[i][0] < max[0])
        {
            for (int x = 0; x < numVecs + 1; ++x)
            {
                max[x] = solutions[i][x];
            }
        }
        avgHeur += solutions[i][0];
        avgTime += executionTime[i];
        csv << i << ", " << (executionTime[i] / 1000.00) << ", " << solutions[i][0] << ", ";
        for (int x = 1; x <= numVecs; ++x)
        {
            csv << labels[solutions[i][x]] << ">";
        }
        csv << labels[0] << "\n";
    }
    avgHeur /= RUN_MANY;
    std::cout << "Average cost: " << avgHeur << std::endl;
    avgTime /= (RUN_MANY * 1000.00);
    std::cout << "Average Time: " << avgTime << std::endl;
    std::cout << "Best Solution:\n\tCost: " << max[0] << "\n\t";
    for (int x = 1; x < numVecs + 1; ++x)
    {
        std::cout << labels[max[x]] << ">";
    }
    std::cout << std::endl;
    csv.close();
    csv.open("ils_iteration_limits.csv");
    csv << "Iterations, Time, Avg Cost, Best Cost \n";
    for (int x = 0; x < MAX_CONSECUTIVE_ITERATIONS; ++x)
    {
        std::cout << (x + 1) << "(" << (avgTimePerIterations[x] / 1000.00) << "): " << avgCostPerIterations[x] << " " << bestCostPerIterations[x] << std::endl;
        csv << (x + 1) << ", " << (avgTimePerIterations[x] / 1000.00) << ", " << avgCostPerIterations[x] << ", " << bestCostPerIterations[x] << "\n";
    }
    csv.close();
#endif
    return 0;
}