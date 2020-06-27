#include "Util.h"
#include "Algorithm.h"

#include <cmath>
#include <algorithm>
#include <libconfig.h++>
#include <assert.h>
#include <float.h>

using namespace libconfig;

using namespace std;

double Util::RandomDouble(double min, double max) {
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}


double Util::gaussRND(double mean, double std) {
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return std * X + mean;
}


bool Util::Load_RC_File(const char *filename) {
    Config cfg;

    cfg.readFile(filename);

    const Setting& root = cfg.getRoot();

    Algorithm::GenerationCount = root["GenerationCount"];
    Algorithm::ChangePeriod = root["ChangePeriod"];
    Algorithm::CrossoverProbability = root["CrossoverProbability"];
    Algorithm::MutationProbability = root["MutationProbability"];
    Algorithm::RandomImmigrantsChangePercentage = root["RandomImmigrantsChangePercentage"];
    Algorithm::PopulationSize = root["PopulationSize"];
    Algorithm::TournamentSize = root["TournamentSize"];
    Algorithm::ExplicitMemorySize = root["ExplicitMemorySize"];
    Algorithm::MemoryUpdateFreq = root["MemoryUpdateFreq"];
    Algorithm::NumberOfPeaks = root["NumberOfPeaks"];
    Algorithm::ShiftLength = root["ShiftLength"];
    Algorithm::SOS_Lambda = root["Lambda"];
    Algorithm::SOS_UseBasisFunction = root["UseBasisFunction"];
    Algorithm::HyperMutationProbability = root["HyperMutationProbability"];
    Algorithm::SOS_ForkingGenerationPeriod = root["ForkingGenerationPeriod"];
    Algorithm::SOS_MinScoutPopulationSizeRelative = root["MinSubPopulationSizeRelative"];
    Algorithm::SOS_MaxScoutPopulationSizeRelative = root["MaxSubPopulationSizeRelative"];
    Algorithm::SOS_MinDiameterRelative = root["MinDiameterRelative"];
    Algorithm::SOS_MaxDiameterRelative = root["MaxDiameterRelative"];
    Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative = root["MinFitnessOfNewForkingPopulationsRelative"];
    Algorithm::SOS_MinFitnessOfExistingForkingPopulationsRelative = root["MinFitnessOfExistingForkingPopulationsRelative"];
    Algorithm::MinCoordinate = root["MinCoord"];
    Algorithm::MaxCoordinate = root["MaxCoord"];
    Algorithm::SOS_MinBasePopulationSizeRelative = root["MinBasePopulationSizeRelative"];
    Algorithm::DimensionSize = root["Dimension"];
    Algorithm::SOS_DiameterReduceFactor = root["SOSDiameterReduceFactor"];
    Algorithm::SOS_ALpha = root["SOSAlpha"];
    Algorithm::RunCount = root["RunCount"];

    return true;
}

std::vector< Individual* > Util::getWorstIndividuals(size_t worstK, std::vector<Individual> &individuals) {
    vector<Individual*> worstIndividuals;

    double worstFitness;
    double fitness;
    Individual* worstIndividualAtCurrentIteration;

    // get worst k individuals of population
    for (int k = 0; k < worstK; k++) {
        worstFitness = DBL_MAX;

        for (vector<Individual>::iterator i = individuals.begin(); i != individuals.end(); i++) {

            vector<Individual*>::iterator found = find( worstIndividuals.begin(), worstIndividuals.end(), &(*i) );

            // already added to worst individuals array so, skip this element
            if (found != worstIndividuals.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness < worstFitness) {
                worstFitness = fitness;
                worstIndividualAtCurrentIteration = &(*i);
            }
        }

        worstIndividuals.push_back( worstIndividualAtCurrentIteration );
    }

    return worstIndividuals;
}

std::vector< Individual* > Util::getBestIndividuals(size_t bestK, std::vector<Individual> &individuals) {
    vector<Individual*> bestIndividuals;

    double bestFitness;
    size_t bestIndex;
    double fitness;
    Individual* tempBest;

    // get worst k individuals of population
    for (int k = 0; k < bestK; k++) {
        bestFitness = -DBL_MAX;

        for (vector<Individual>::iterator i = individuals.begin();
             i != individuals.end(); i++) {

            vector<Individual*>::iterator found = find( bestIndividuals.begin(), bestIndividuals.end(), &(*i) );

            // already added to worst individuals array so, skip this element
            if (found != bestIndividuals.end())
                continue;

            fitness = (*i).getFitness();

            if (fitness > bestFitness) {
                bestFitness = fitness;
                tempBest = &(*i);
            }
        }

        bestIndividuals.push_back(tempBest);
    }

    return bestIndividuals;
}

double Util::GetEuclideanDistance(const Individual &from, const Individual &to) {
    double* fromValues = from.getValues();
    double* toValues = to.getValues();

    assert( fromValues != nullptr && toValues != nullptr );

    double d = 0;
    for (int i = 0; i < Individual::getDimension(); i++) {
        d += pow(fromValues[i]-toValues[i], 2.0);
    }

    return sqrt( d );
}

void Util::DisplayPercentage(size_t step, size_t total, PercentageType pt) {
    if (pt == PT_Arrow) {
        cout << "[" << flush;
        for (int i = 0; i < total; i++) {
            if (i == step)
                cout << ">" << flush;
            else if (i < step)
                cout << "=" << flush;
            else
                cout << "." << flush;
        }
        cout << "]" << endl;
    }
    else if (pt == PT_Dot) {
        cout << "." << flush;

        if ( step == total-1 ) {
            cout << endl;
        }
    }
}

void Util::BubbleSort(double *arr, int *sortedIndices, size_t size) {
    int pass,j,hold,in;
    for(pass=1;pass<=size-1;pass++){
        for(j=0;j<=size-2;j++){
            if(arr[j]<arr[j+1]){
                hold = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = hold;
                in = sortedIndices[j];
                sortedIndices[j]=sortedIndices[j+1];
                sortedIndices[j+1]=in;
            }
        }
    }
}
