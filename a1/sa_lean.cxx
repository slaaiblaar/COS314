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
#ifdef MAX_ITERATIONS
const int MAX_CONSECUTIVE_ITERATIONS = MAX_ITERATIONS;
#else
const int MAX_CONSECUTIVE_ITERATIONS = 50;
#endif
const int NUM_NEIGHBOURS = 30;
const int temp0 = 1000;
int heuristic(int *sol, int len)
{
  int totalVal = 0;
  for (int x = 0; x < len - 1; ++x)
  {
    totalVal += distances[sol[x]][sol[x + 1]];
  }
  totalVal += distances[sol[len - 1]][sol[0]];
  return totalVal;
}
bool swapped[numVecs][numVecs];
int mutations[NUM_NEIGHBOURS][2];
void howdyNeighbours(int *sol, int neighbours[][numVecs])
{
  int v1, v2;
  // ensures neighbours aren't duplicates
  bool duplicates;
  for (int n = 0; n < NUM_NEIGHBOURS; ++n)
  {
    swapped[mutations[n][0]][mutations[n][1]] = false;
    swapped[mutations[n][1]][mutations[n][0]] = false;
  }
  for (int n = 0; n < NUM_NEIGHBOURS; ++n)
  {
    do
    {
      v1 = (rand() % (numVecs - 1)) + 1;
      v2 = (rand() % (numVecs - 1)) + 1;
    } while (v1 == v2 && (swapped[v1][v2] || swapped[v1][v2]));
    for (int v = 0; v < numVecs; ++v)
    {
      neighbours[n][v] = sol[v];
    }
    neighbours[n][v1] = sol[v2];
    neighbours[n][v2] = sol[v1];
    swapped[v1][v2] = true;
    swapped[v2][v1] = true;
    mutations[n][0] = v1;
    mutations[n][1] = v2;
  }
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
  // Don't want reallocation of memory to factor into execution time
  // when running multiple iterations for data gathering
  double temp;
  int consecutiveIterations;
  int currSol[numVecs];
  int bestSol[numVecs];
  int neighbours[NUM_NEIGHBOURS][numVecs];
  int currHeur;
  int bestHeur;
  int neighbourHeur;
  int bestNeighbour;
  double random;
  std::chrono::_V2::system_clock::time_point time1; // disgusting
  std::chrono::_V2::system_clock::time_point time2; // disgusting
  int maxConsecutiveIterations = MAX_CONSECUTIVE_ITERATIONS;
#ifdef RUN_MANY
  int iterations = 0;
  int solutions[RUN_MANY][numVecs + 2];
  int randomPadding[RUN_MANY][numVecs + 2];
  double avgNeighbourHeuristics[RUN_MANY][MAX_CONSECUTIVE_ITERATIONS];
  int moreRandomPadding[RUN_MANY][MAX_CONSECUTIVE_ITERATIONS];
  long int executionTime[RUN_MANY];
  long int evenMoreRandomPadding[RUN_MANY];
  int totalDiffCurrBest[RUN_MANY] = {0}; // Sum of (curHeur - bestHeur) for all iterations
  maxConsecutiveIterations = 1;
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
    avgCostPerIterations[maxConsecutiveIterations - 1] /= RUN_MANY * 1.0;
    avgTimePerIterations[maxConsecutiveIterations - 1] /= RUN_MANY * 1.0;
    std::cout << maxConsecutiveIterations << ": " << avgCostPerIterations[maxConsecutiveIterations - 1] << " cost, " << (avgTimePerIterations[maxConsecutiveIterations - 1] / 1000.0) << "ms avg\n";
    std::cout << "Max iterations: " << maxConsecutiveIterations << std::endl;
    ++maxConsecutiveIterations;
    iterations = 0;
  }
  if (maxConsecutiveIterations > MAX_CONSECUTIVE_ITERATIONS)
    goto end_loop;
#endif
  time1 = std::chrono::high_resolution_clock::now();
  temp = temp0;
  consecutiveIterations = 0; // t

#if SEARCH_SPACE == LARGE_SPACE
  calculateVecDists();
#endif
  std::fill(currSol, currSol + numVecs, 0);

  currSol[0] = 0;
  // Construct an initial solution by generating a random number, n, in range [0, number of unvisited nodes]
  // and then adding the nth unvisited node to the solution
  // only works because the solution space is a bi-directional, connected graph
  // otherwise a constructive method would need to be used
  generateInitialSol(currSol);

  std::cout << "Initial Solution: ";
  for (int x = 0; x < numVecs; ++x)
  {
    bestSol[x] = currSol[x];
    std::cout << labels[currSol[x]] << "->";
  }
  std::cout << labels[0] << std::endl;
  std::cout << "Heuristic Value: " << heuristic(currSol, numVecs) << "\n";
  currHeur = heuristic(currSol, numVecs);
  bestHeur = heuristic(bestSol, numVecs);
  while (consecutiveIterations < maxConsecutiveIterations)
  {
    ++consecutiveIterations;

    howdyNeighbours(currSol, neighbours);
    bestNeighbour = 0;
    neighbourHeur = heuristic(neighbours[0], numVecs);
    for (int n = 1; n < NUM_NEIGHBOURS; ++n)
    {
      if (heuristic(neighbours[n], numVecs) < neighbourHeur)
      {
        bestNeighbour = n;
        neighbourHeur = heuristic(neighbours[n], numVecs);
      }
    }
    random = drand48();
    std::cout << "\tBest neighbour\n\t\t";
    for (int x = 0; x < numVecs; ++x)
    {
      std::cout << labels[neighbours[bestNeighbour][x]] << "->";
    }
    std::cout << labels[neighbours[bestNeighbour][0]] << std::endl;
    if (random < exp((currHeur - neighbourHeur) / temp))
    {
      std::cout << "\tSwitching currSol to neighbour\n";
      for (int x = 0; x < numVecs; ++x)
      {
        currSol[x] = neighbours[bestNeighbour][x];
      }
      currHeur = neighbourHeur;
      std::cout << "\t\tBoltzmann factor: " << exp((currHeur - neighbourHeur) / temp) << std::endl;
      std::cout << "\t\tRandom[0,1): " << random << std::endl;
    }
    if (currHeur < bestHeur)
    {
      std::cout << "\tNew best sol: \n\t\t";
      for (int x = 0; x < numVecs; ++x)
      {
        std::cout << labels[bestSol[x] = currSol[x]] << "->";
      }
      std::cout << labels[0] << std::endl;
      std::cout << "\t\tHeuristic Value: " << (bestHeur = currHeur) << std::endl;
    }
    temp = temp0 / log(consecutiveIterations + 1);
  }
  time2 = std::chrono::high_resolution_clock::now();
#ifndef RUN_MANY // Outputs the best solution if the algorithm was only being run once
  std::cout << std::endl
            << (std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1000.00) << ":Optimal solution: \n";
  for (int x = 0; x < numVecs; ++x)
  {
    std::cout << labels[bestSol[x]] << "->";
  }
  std::cout << std::endl
            << labels[0]
            << "\nCost: " << heuristic(bestSol, numVecs) << std::endl;
#endif
#ifdef RUN_MANY // if algorithm was run multiple times, write statistics to files and display data
  executionTime[iterations] = std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count();
  solutions[iterations][0] = bestHeur;
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
  csv.open("sa_data.csv");
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
      csv << labels[solutions[i][x]] << "->";
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
  csv.open("sa_iteration_limits.csv");
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