#ifndef EVOLUTIONARYLIB_ALGORITHM_H
#define EVOLUTIONARYLIB_ALGORITHM_H


#include <glob.h>
#include "Population.h"
#include "ScoutPopulation.h"
#include "ParentPopulation.h"
#include "Statistics.h"

class Algorithm {
public:
    // Crossover type
    enum CX_TYPE {
    SBX, HillClimbing
    };

    static unsigned int PopulationSize;
    static unsigned int GenerationCount;
    static unsigned int ChangePeriod;
    static double CrossoverProbability;
    static double MutationProbability;
    static unsigned int TournamentSize;
    static bool Elitism;
    static unsigned int RunCount;
    static Statistics statistics;

    // Random Immigrants parameters
    static double RandomImmigrantsChangePercentage;

    // Memory Search parameters
    static unsigned int ExplicitMemorySize;
    static unsigned int MemoryUpdateFreq;

    // Hyper mutation parameters
    static double HyperMutationProbability;

    // Moving Peaks Parameters
    static unsigned int NumberOfPeaks;
    static double MinCoordinate;
    static double MaxCoordinate;
    static unsigned int DimensionSize;
    static double ShiftLength;


    // SOS parameters
    static double SOS_Lambda;
    static bool   SOS_UseBasisFunction;
    static unsigned int SOS_ForkingGenerationPeriod;
    static double SOS_MinScoutPopulationSizeRelative;
    static double SOS_MaxScoutPopulationSizeRelative;
    static double SOS_MinBasePopulationSizeRelative;
    static double SOS_MinDiameterRelative;
    static double SOS_MaxDiameterRelative;
    static double SOS_MinFitnessOfNewForkingPopulationsRelative;
    static double SOS_MinFitnessOfExistingForkingPopulationsRelative;
    static double SOS_DiameterReduceFactor;
    static double SOS_ALpha;
    static CX_TYPE CX_Type;

    static void InitMovPeaks();
    static void ReleaseMovPeaks();


    static void SelectParents( Population* population, int* momIndex, int* dadIndex, size_t tournamentSize);
    // static void SelectParentsRankBased( Population* population, int* momIndex, int* dadIndex, size_t tournamentSize);
    // static void AssignRanks( Population* population );
    static std::vector<Individual> SimulatedBinaryCrossover( const Individual& mom, const Individual& dad );
    static std::vector<Individual> HillClimbingCrossover( Individual& mom, Individual& dad );
    static void GenerationalGA( Population* population );
    static void GenerationalGA( Population* population, size_t suggestedPopulationSize, CX_TYPE cxType=SBX);

    static void GenerationalGARandomImmigrants( Population* population );
    static void ApplyElitism( Population* population, const Individual& best );

    // simple ga
    static void RunSimpleGA();

    // hyper mutation
    static void RunHyperMutation();

    // random immigranats
    static void RunRandomImmigrants();


    // memory search
    static void RunMemorySearch();
    static void UpdateMemoryPopulation(Population &memoryPopulation, vector<Individual> &explicitMemory);
    static void UpdateExplicitMemory( Population& memoryPopulation,
                                      Population& searchPopulation, std::vector<Individual>& explicitMemory );
    static Individual* MinDistIndividual( std::vector<Individual>& explicitMemory, Individual* bestIndividual );

    // self-organizing scouts
    static void RunSelfOrganizingScouts();
    static void RunSelfOrganizingScoutsRandomImmigrants(bool fix);
    static void SOSGenerationalGA( ParentPopulation* parentPopulation );
    static void SOSGenerationalGARandomImmigrants( ParentPopulation* parentPopulation );
    static void SOSAdjustPopulationSizes(ParentPopulation *parentPopulation);
    static void SOSAdjustPopulationSizesRandomImmigrants(ParentPopulation *parentPopulation);
    static void SOSAdjustSearchSpace( ParentPopulation* parentPopulation );
    static void SOSDeleteNonFitScouts( ParentPopulation* parentPopulation );
    static double SOSGetRelativeDynamism(Population *population, double totalDynamism);
    static double SOSGetRelativeFitness(Population *population, double overallMinFitness, double totalFitness);
    static double SOSGetRelativeQuality(Population *population,
                                        double beta, double totalDynamism, double overallMinFitness,
                                        double totalFitness);
    static double SOSGetQuality(Population *population,
                                double beta, double totalDynamism, double overallMinFitness,
                                double totalFitness);
    static size_t SOSSuggestedSize(Population *population, double quality);
    static ScoutPopulation SOSFork(ParentPopulation *parentPopulation);
    static ScoutPopulation SOSForkSupposedToBe(ParentPopulation *parentPopulation);
    //static ScoutPopulation SOSForkDemet(ParentPopulation *parentPopulation);
    static bool SOSMergeScoutPopulations( ParentPopulation *parentPopulation );
    static void MergeScoutPopulations( ParentPopulation* parentPopulation, ScoutPopulation* better, ScoutPopulation* worse );

    static void SOSReduceScoutDiameters( ParentPopulation* parentPopulation );

    static void SOSAddRandomImmigrants( Population* population, size_t size );
    static void SOSRemoveWorstIndividuals( Population* population, size_t size );


    // moving center sos
    static void GenerationalGAMovingCenter( Population* population );
    static void RunSelfOrganizingScoutsRandomImmigrantsMovingCenter();
    static void SOSGenerationalGARandomImmigrantsMovingCenter( ParentPopulation* parentPopulation );

};


#endif //EVOLUTIONARYLIB_ALGORITHM_H
