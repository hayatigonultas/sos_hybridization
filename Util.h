#ifndef EVOLUTIONARYLIB_UTIL_H
#define EVOLUTIONARYLIB_UTIL_H


#include <stddef.h>
#include <iostream>
#include <vector>
#include <iomanip>

class Individual;
class Population;

using std::cout;
using std::flush;
using std::setw;

class Util {
public:
    enum PercentageType { PT_Arrow, PT_Dot, PT_Percent };
    static double RandomDouble(double min, double max);

    static double gaussRND(double mean, double std);

    static bool Load_RC_File(const char* filename);

    template <typename TYPE>
    static void printArray(TYPE* arr, size_t size) {
        for (int i = 0; i < size; i++) {
            std::cout << arr[i] << std::endl;
        }
    }

    static std::vector< Individual* > getWorstIndividuals(size_t worstK, std::vector<Individual>& individuals);
    static std::vector< Individual* > getBestIndividuals(size_t bestK, std::vector<Individual>& individuals);

    static double GetEuclideanDistance( const Individual& from, const Individual& to );

    static void DisplayPercentage(size_t step, size_t total, PercentageType pt = PT_Arrow);

    static void BubbleSort( double* arr, int* sortedIndices, size_t size);

    static inline void LoadBar(unsigned int x, unsigned int n, unsigned int w = 50)
    {
        if ( (x != n) && (x % (n/100+1) != 0) ) return;

        float ratio  =  x/(float)n;
        int   c      =  ratio * w;

        cout << setw(3) << (int)(ratio*100) << "% [";
        for (int x=0; x<c; x++) cout << "=";
        for (int x=c; x<w; x++) cout << " ";
        cout << "]\r" << flush;
    }
};


#endif //EVOLUTIONARYLIB_UTIL_H
