#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <time.h>
const int space = 50;
int main()
{
    srand(time(NULL));
    int numExperiments = 20000;

    int freq[space] = {0};

    for (int x = 0; x < numExperiments; ++x)
    {
        ++freq[rand() % space];
    }

    for (int x = 0; x < space; ++x)
    {
        std::cout << x << ":\t" << freq[x] << "\t";
        int tally = 0;
        do
        {
            std::cout << "|";
            tally += 5;
        } while (tally < freq[x]);
        std::cout << std::endl;
    }

    return 0;
}