/*
 * Hayati Gonultas
 * All rights reserved
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "Algorithm.h"
#include "movpeaks.h"
#include "Util.h"

using namespace std;


unsigned int Algorithm::PopulationSize;
unsigned int Algorithm::GenerationCount;
unsigned int Algorithm::ChangePeriod;
double Algorithm::CrossoverProbability;
double Algorithm::MutationProbability;
unsigned int Algorithm::TournamentSize;
bool Algorithm::Elitism = true;
unsigned int Algorithm::RunCount;
Statistics Algorithm::statistics;

// Random Immigrants parameters
double Algorithm::RandomImmigrantsChangePercentage;

// Memory Search parameters
unsigned int Algorithm::ExplicitMemorySize;
unsigned int Algorithm::MemoryUpdateFreq;

// Hyper mutation parameters
double Algorithm::HyperMutationProbability;

// Moving Peaks Parameters
unsigned int Algorithm::NumberOfPeaks;
double Algorithm::MinCoordinate;
double Algorithm::MaxCoordinate;
unsigned int Algorithm::DimensionSize;
double Algorithm::ShiftLength;

// SOS parameters
double Algorithm::SOS_Lambda;
bool   Algorithm::SOS_UseBasisFunction;
unsigned int Algorithm::SOS_ForkingGenerationPeriod;
double Algorithm::SOS_MinScoutPopulationSizeRelative;
double Algorithm::SOS_MaxScoutPopulationSizeRelative;
double Algorithm::SOS_MinBasePopulationSizeRelative;
double Algorithm::SOS_MinDiameterRelative;
double Algorithm::SOS_MaxDiameterRelative;
double Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative;
double Algorithm::SOS_MinFitnessOfExistingForkingPopulationsRelative;
double Algorithm::SOS_DiameterReduceFactor;
double Algorithm::SOS_ALpha;
Algorithm::CX_TYPE Algorithm::CX_Type = Algorithm::SBX;

extern int number_of_peaks;
extern double vlength;
extern double lambda;
extern int use_basis_function;
extern double mincoordinate;
extern double maxcoordinate;

void Algorithm::InitMovPeaks() {
    number_of_peaks = Algorithm::NumberOfPeaks;
    vlength = Algorithm::ShiftLength;
    lambda = Algorithm::SOS_Lambda;
    use_basis_function = Algorithm::SOS_UseBasisFunction;
    mincoordinate = Algorithm::MinCoordinate;
    maxcoordinate = Algorithm::MaxCoordinate;

    init_peaks();
}

//void Algorithm::SelectParentsRankBased(Population *population, int *momIndex, int *dadIndex, size_t tournamentSize) {
//    size_t size = std::min( population->getPopulationSize(), tournamentSize ); // population size may be less than tournament size
//
//    vector<int> indices(size, -1);
//
//    vector<Individual>& individuals = population->getIndividuals();
//
//    unsigned int newGeneratedValue;
//
//    // get indices
//    for (int i = 0; i < size; ) {
//        newGeneratedValue = (unsigned int) (rand() % population->getPopulationSize());
//
//        vector<int>::iterator it = find(indices.begin(), indices.end(), newGeneratedValue);
//        if (it != indices.end()) { // we already added this index to vector
//            continue;
//        }
//
//        indices[i] = newGeneratedValue;
//        i++;
//    }
//
//    // get momIndex
//    double bestFitness = -DBL_MAX;
//    for (int i = 0; i < size; i++) {
//        double f = individuals[indices[i]].getFitness();
//
//        if (f > bestFitness) {
//            bestFitness = f;
//            *momIndex = indices[i];
//        }
//    }
//
//    // get dadIndex
//    bestFitness = -DBL_MAX;
//    for (int i = 0; i < size; i++) {
//        if (indices[i] == *momIndex) { // this is mom index, so skip
//            continue;
//        }
//
//        double f = individuals[indices[i]].getFitness();
//
//        if (f > bestFitness) {
//            bestFitness = f;
//            *dadIndex = indices[i];
//        }
//    }
//}

void Algorithm::SelectParents(Population* population, int *momIndex, int *dadIndex, size_t tournamentSize) {
    size_t size = std::min( population->getPopulationSize(), tournamentSize ); // population size may be less than tournament size

    vector<int> indices(size, -1);

    vector<Individual>& individuals = population->getIndividuals();

    unsigned int newGeneratedValue;

    // get indices
    for (int i = 0; i < size; ) {
        newGeneratedValue = (unsigned int) (rand() % population->getPopulationSize());

        vector<int>::iterator it = find(indices.begin(), indices.end(), newGeneratedValue);
        if (it != indices.end()) { // we already added this index to vector
            continue;
        }

        indices[i] = newGeneratedValue;
        i++;
    }

    // get momIndex
    double bestFitness = -DBL_MAX;
    for (int i = 0; i < size; i++) {
        double f = individuals[indices[i]].getFitness();

        if (f > bestFitness) {
            bestFitness = f;
            *momIndex = indices[i];
        }
    }

    // get dadIndex
    bestFitness = -DBL_MAX;
    for (int i = 0; i < size; i++) {
        if (indices[i] == *momIndex) { // this is mom index, so skip
            continue;
        }

        double f = individuals[indices[i]].getFitness();

        if (f > bestFitness) {
            bestFitness = f;
            *dadIndex = indices[i];
        }
    }
}



void Algorithm::ReleaseMovPeaks() {
    free_peaks();
}


void Algorithm::GenerationalGAMovingCenter(Population *population) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    assert( population->getPopulationSize() >= 4 );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < population->getPopulationSize() && (i+1) < population->getPopulationSize(); i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());
            offspring[1].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividualMovingCenter(offspring[0])
                || !population->isValidIndividualMovingCenter(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( population->getPopulationSize() % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividualMovingCenter(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}

void Algorithm::GenerationalGA(Population* population) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    assert( population->getPopulationSize() >= 4 );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < population->getPopulationSize() && (i+1) < population->getPopulationSize(); i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());
            offspring[1].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]) || !population->isValidIndividual(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( population->getPopulationSize() % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}

void Algorithm::GenerationalGARandomImmigrants(Population *population) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < population->getPopulationSize() && (i+1) < population->getPopulationSize(); i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());
            offspring[1].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]) || !population->isValidIndividual(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( population->getPopulationSize() % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );

            offspring[0].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}


void Algorithm::RunSimpleGA() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) {
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();


            change_peaks();
            population.generatePopulationRandom();

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
        }

        GenerationalGA( &population );

        statistics.addStat( population.best().getFitness() );
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


void Algorithm::ApplyElitism(Population* population, const Individual &best) {
    if (Algorithm::Elitism) {
        int index;
        population->worst(&index);
        population->setIndividual( best, index );
    }
}

void Algorithm::RunRandomImmigrants() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    size_t randomImmigrantsSize = (size_t) (Algorithm::PopulationSize * Algorithm::RandomImmigrantsChangePercentage);

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) {
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            change_peaks();
            population.updateFitnesses();

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
        }

        vector<Individual>& individuals = population.getIndividuals();

        vector<Individual*> worstIndividuals =
                Util::getWorstIndividuals(randomImmigrantsSize, individuals);

        for(vector<Individual*>::iterator it = worstIndividuals.begin(); it != worstIndividuals.end(); it++) {
            population.initIndividualRandom( *it );
        }

        GenerationalGA( &population );

        statistics.addStat( population.best().getFitness() );
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::RunMemorySearch() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population searchPopulation( 45 );
    Population memoryPopulation( 45 );

    vector<Individual> explicitMemory;

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) { // env. change
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            // merge explicit mem and memory pop select best n
            change_peaks();

            // re-init search population
            searchPopulation.generatePopulationRandom();

            // refresh fitnesses of memory individuals
            memoryPopulation.updateFitnesses();

            statistics.addBestErrorAtChangeStat( memoryPopulation.best().getFitness() );

            // refresh fitnesses of explicit memory
            for (int j = 0; j < explicitMemory.size(); j++) {
                explicitMemory[j].updateFitness();
            }

            // update memory population
            Algorithm::UpdateMemoryPopulation( memoryPopulation, explicitMemory );
        }

        if (i > 0 && (i % Algorithm::MemoryUpdateFreq == 0)) {
            // update explicit memory
            Algorithm::UpdateExplicitMemory( memoryPopulation, searchPopulation, explicitMemory );
        }

        // generational part of memory population
        GenerationalGA( &memoryPopulation );

        // generational part of search population
        GenerationalGA( &searchPopulation );

        statistics.addStat( memoryPopulation.best().getFitness() );
        statistics.addAbsoluteRecoveryRateStat( memoryPopulation.best().getFitness() );
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::UpdateMemoryPopulation(Population &memoryPopulation, vector<Individual> &explicitMemory) {
    vector<Individual> mergedIndividuals = memoryPopulation.getIndividuals();
    mergedIndividuals.insert( mergedIndividuals.end(), explicitMemory.begin(), explicitMemory.end() );

    vector<Individual*> bestIndividuals = Util::getBestIndividuals( memoryPopulation.getPopulationSize(), mergedIndividuals );

    // copy best individuals to memory population
    for (int i = 0; i < memoryPopulation.getPopulationSize(); i++) {
        memoryPopulation.setIndividual( *(bestIndividuals[i]), i );
    }
}

void Algorithm::UpdateExplicitMemory(Population &memoryPopulation, Population &searchPopulation,
                                     std::vector<Individual> &explicitMemory) {
    Individual bestIndividualSearch = searchPopulation.best();
    Individual bestIndividualMemory = memoryPopulation.best();

    Individual* bestIndividual;

    if (bestIndividualMemory.getFitness() > bestIndividualSearch.getFitness())
        bestIndividual = &bestIndividualMemory;
    else
        bestIndividual = &bestIndividualSearch;


    if (explicitMemory.size() < Algorithm::ExplicitMemorySize) { // explicit memory with free slots
        explicitMemory.push_back( *bestIndividual );
    }
    else {
        // memory is full find mindist candidate to replace
        Individual* mindistIndividual = MinDistIndividual( explicitMemory, bestIndividual );

        // update explicit memory
        if (mindistIndividual->getFitness() < bestIndividual->getFitness()) {
            *mindistIndividual = *bestIndividual;
        }
    }
}

Individual *Algorithm::MinDistIndividual(std::vector<Individual> &explicitMemory, Individual *bestIndividual) {
    double mindist = DBL_MAX;
    Individual* mindistIndividual = nullptr;
    double distance;

    for (int i = 0; i < explicitMemory.size(); i++) {
        distance = bestIndividual->getEuclideanDistanceToIndividual(explicitMemory[i]);
        if ( distance < mindist ) {
            mindist = distance;
            mindistIndividual = &(explicitMemory[i]);
        }
    }

    return  mindistIndividual;
}

void Algorithm::RunSelfOrganizingScouts() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    ParentPopulation parentPopulation( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        SOSGenerationalGA( &parentPopulation );

        SOSAdjustSearchSpace( &parentPopulation );

        SOSDeleteNonFitScouts( &parentPopulation );

        SOSAdjustPopulationSizes( &parentPopulation );

        ScoutPopulation sp = SOSFork( &parentPopulation );

        if (sp.getIndividuals().size() > 0) { // forked a child
            parentPopulation.addScoutPopulation( sp );
            parentPopulation.generatePopulationRandom();
        }

        while (SOSMergeScoutPopulations( &parentPopulation )) {
        }

        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) { // env. change
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            // merge explicit mem and memory pop select best n
            change_peaks();

            parentPopulation.updateAllFitnesses();

            SOSAdjustSearchSpace( &parentPopulation );

            Algorithm::statistics.addBestErrorAtChangeStat( parentPopulation.overallBestFitness() );
        }

        assert( parentPopulation.getPopulationSize() >= 8 );

        statistics.addStat( parentPopulation.overallBestFitness() );
        statistics.addAbsoluteRecoveryRateStat( parentPopulation.overallBestFitness() );

        SOSReduceScoutDiameters( &parentPopulation );
//        vector<ScoutPopulation>& scouts = parentPopulation.getScoutPopulations();
//        for(int i = 0; i < scouts.size(); i++) {
//                assert( scouts[i].getPopulationSize() >= 4 &&  scouts[i].getPopulationSize() <= 10);
//        }
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


void Algorithm::SOSGenerationalGA(ParentPopulation *parentPopulation) {
    // generational ga steps for
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); ++it) {
        if (CX_Type == HillClimbing) {
            GenerationalGA( &(*it), (*it).getSuggestedSize(), HillClimbing );
        }
        else {
            GenerationalGA( &(*it), (*it).getSuggestedSize(), SBX );
        }
    }

    GenerationalGA( parentPopulation, parentPopulation->getSuggestedSize(), SBX );
}

void Algorithm::SOSGenerationalGARandomImmigrants(ParentPopulation *parentPopulation) {
    // generational ga steps for
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); ++it) {
        GenerationalGA( &(*it)  );
    }

    GenerationalGA( parentPopulation );
}

void Algorithm::SOSAdjustPopulationSizesRandomImmigrants(ParentPopulation *parentPopulation) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);


    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();


    double obf = parentPopulation->overallBestFitness();

    // parent population (1) + number of child populations
    size_t numberOfPopulations = scouts.size() + 1;

    // get fit populations count and minimum fitness of all individuals (populations whose fitnesses are incrementing)
    size_t fitPopulationCount = 0;
    double minFitness = DBL_MAX;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {

        // get min fitness
        if ( minFitness > (*it).getFitness() ) {
            minFitness = (*it).getFitness();
        }

        // is incrementing?
        if ((*it).getFitness() > (*it).getPreviousBest().getFitness() ) {
            // will remove this population
            fitPopulationCount++;
        }
    }

    // get min fitness compare base population
    if ( minFitness > parentPopulation->getFitness() ) {
        minFitness = parentPopulation->getFitness();
    }

    // is incrementing for parent?
    // is incrementing?
    if (parentPopulation->getFitness() > parentPopulation->getPreviousBest().getFitness() ) {
        // will remove this population
        fitPopulationCount++;
    }

    double beta = (double) fitPopulationCount / numberOfPopulations;

    // calculate sumFitness
    double sumFitness = 0;
    double sumDynamism = 0;
    sumFitness += parentPopulation->getFitness() - minFitness;
    sumDynamism += parentPopulation->getDynamism();
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumFitness += (*it).getFitness() - minFitness;
        sumDynamism += (*it).getDynamism();
    }

    // arrange individuals between populations
    // first check out populations that is offered individual count > max individual count
    vector<ScoutPopulation*> minScoutSizeScouts;
    vector<ScoutPopulation*> maxScoutSizeScouts;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        double quality = Algorithm::SOSGetQuality(&(*it), beta, sumDynamism, minFitness, sumFitness);
        (*it).setQuality( quality );
    }

    double parentRelDynamism = SOSGetRelativeDynamism(parentPopulation, sumDynamism);
    double parentRelFitness = SOSGetRelativeFitness(parentPopulation, minFitness, sumFitness);
    double parentQuality = SOSGetRelativeQuality(parentPopulation, beta, sumDynamism, minFitness, sumFitness);

    double sumQuality = parentQuality;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumQuality += (*it).getQuality();
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setQuality( (*it).getQuality()/sumQuality );
    }

    parentQuality /= sumQuality;

    size_t totalScoutSizes = 0;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setSuggestedSize( SOSSuggestedSize( &(*it), (*it).getQuality() ) );
        totalScoutSizes += (*it).getSuggestedSize();
    }

    parentPopulation->setSuggestedSize( Algorithm::PopulationSize - totalScoutSizes );

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() > maxScoutSize ) {
            parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize() + (*it).getSuggestedSize()-maxScoutSize );
            (*it).setSuggestedSize( maxScoutSize );
        }
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() < minScoutSize ) {
            size_t maximumSizeScout = 0;
            ScoutPopulation* maximumScout = nullptr;
            for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
                if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                    maximumSizeScout = (*itMax).getSuggestedSize();
                    maximumScout = &(*itMax);
                }
            }

            if (maximumScout != nullptr && maximumScout->getSuggestedSize() > parentPopulation->getSuggestedSize()) {
                maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }
            else {
                parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }

            (*it).setSuggestedSize( minScoutSize );
        }
    }

    if (parentPopulation->getSuggestedSize() < minBaseSize) {
        size_t maximumSizeScout = 0;
        ScoutPopulation* maximumScout = nullptr;
        for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
            if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                maximumSizeScout = (*itMax).getSuggestedSize();
                maximumScout = &(*itMax);
            }
        }

        maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minBaseSize-parentPopulation->getSuggestedSize()) );
        parentPopulation->setSuggestedSize( minBaseSize );
    }


    size_t totalIndividualsOfScouts = parentPopulation->getSuggestedSize();
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        totalIndividualsOfScouts += (*it).getSuggestedSize();
    }

    cout << "total: " << totalIndividualsOfScouts << endl;

    assert( totalIndividualsOfScouts == Algorithm::PopulationSize );


    // random immmigrants part

    // parent population
    if (parentPopulation->getPopulationSize() < parentPopulation->getSuggestedSize()) {
        SOSAddRandomImmigrants( parentPopulation, parentPopulation->getSuggestedSize()-parentPopulation->getPopulationSize() );
    }
    else if (parentPopulation->getPopulationSize() > parentPopulation->getSuggestedSize()) {
        SOSRemoveWorstIndividuals( parentPopulation, parentPopulation->getPopulationSize()-parentPopulation->getSuggestedSize() );
    }

    // scouts
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ((*it).getPopulationSize() < (*it).getSuggestedSize()) {
            SOSAddRandomImmigrants( &(*it), (*it).getSuggestedSize()-(*it).getPopulationSize() );
        }
        else if ((*it).getPopulationSize() > (*it).getSuggestedSize()) {
            SOSRemoveWorstIndividuals( &(*it), (*it).getPopulationSize()-(*it).getSuggestedSize() );
        }
    }
}

void Algorithm::SOSAdjustPopulationSizes(ParentPopulation *parentPopulation) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);


    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();


    double obf = parentPopulation->overallBestFitness();

    // parent population (1) + number of child populations
    size_t numberOfPopulations = scouts.size() + 1;

    // get fit populations count and minimum fitness of all individuals (populations whose fitnesses are incrementing)
    size_t fitPopulationCount = 0;
    double minFitness = DBL_MAX;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {

        // get min fitness
        if ( minFitness > (*it).getFitness() ) {
            minFitness = (*it).getFitness();
        }

        // is incrementing?
        if ((*it).getFitness() > (*it).getPreviousBest().getFitness() ) {
            // will remove this population
            fitPopulationCount++;
        }
    }

    // get min fitness compare base population
    if ( minFitness > parentPopulation->getFitness() ) {
        minFitness = parentPopulation->getFitness();
    }

    // is incrementing for parent?
    // is incrementing?
    if (parentPopulation->getFitness() > parentPopulation->getPreviousBest().getFitness() ) {
        // will remove this population
        fitPopulationCount++;
    }

    double beta = (double) fitPopulationCount / numberOfPopulations;

    // calculate sumFitness
    double sumFitness = 0;
    double sumDynamism = 0;
    sumFitness += parentPopulation->getFitness() - minFitness;
    sumDynamism += parentPopulation->getDynamism();
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumFitness += (*it).getFitness() - minFitness;
        sumDynamism += (*it).getDynamism();
    }

    // arrange individuals between populations
    // first check out populations that is offered individual count > max individual count
    vector<ScoutPopulation*> minScoutSizeScouts;
    vector<ScoutPopulation*> maxScoutSizeScouts;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        double quality = Algorithm::SOSGetQuality(&(*it), beta, sumDynamism, minFitness, sumFitness);
        (*it).setQuality( quality );
    }

    double parentRelDynamism = SOSGetRelativeDynamism(parentPopulation, sumDynamism);
    double parentRelFitness = SOSGetRelativeFitness(parentPopulation, minFitness, sumFitness);
    double parentQuality = SOSGetRelativeQuality(parentPopulation, beta, sumDynamism, minFitness, sumFitness);

    double sumQuality = parentQuality;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        sumQuality += (*it).getQuality();
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setQuality( (*it).getQuality()/sumQuality );
    }

    parentQuality /= sumQuality;

    size_t totalScoutSizes = 0;
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        (*it).setSuggestedSize( SOSSuggestedSize( &(*it), (*it).getQuality() ) );
        totalScoutSizes += (*it).getSuggestedSize();
    }

    parentPopulation->setSuggestedSize( Algorithm::PopulationSize - totalScoutSizes );

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() > maxScoutSize ) {
            parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize() + (*it).getSuggestedSize()-maxScoutSize );
            (*it).setSuggestedSize( maxScoutSize );
        }
    }

    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        if ( (*it).getSuggestedSize() < minScoutSize ) {
            size_t maximumSizeScout = 0;
            ScoutPopulation* maximumScout = nullptr;
            for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
                if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                    maximumSizeScout = (*itMax).getSuggestedSize();
                    maximumScout = &(*itMax);
                }
            }

            if (maximumScout != nullptr && maximumScout->getSuggestedSize() > parentPopulation->getSuggestedSize()) {
                maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }
            else {
                parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-(minScoutSize-(*it).getSuggestedSize()) );
            }

            (*it).setSuggestedSize( minScoutSize );
        }
    }

    if (parentPopulation->getSuggestedSize() < minBaseSize) {
        size_t maximumSizeScout = 0;
        ScoutPopulation* maximumScout = nullptr;
        for ( std::vector<ScoutPopulation>::iterator itMax = scouts.begin(); itMax != scouts.end(); itMax++) {
            if ( (*itMax).getSuggestedSize() > maximumSizeScout ) {
                maximumSizeScout = (*itMax).getSuggestedSize();
                maximumScout = &(*itMax);
            }
        }

        maximumScout->setSuggestedSize( maximumScout->getSuggestedSize()-(minBaseSize-parentPopulation->getSuggestedSize()) );
        parentPopulation->setSuggestedSize( minBaseSize );
    }


    size_t totalIndividualsOfScouts = parentPopulation->getSuggestedSize();
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        totalIndividualsOfScouts += (*it).getSuggestedSize();
    }

//    cout << "total: " << totalIndividualsOfScouts << endl;

    assert( totalIndividualsOfScouts == Algorithm::PopulationSize );
}

void Algorithm::SOSDeleteNonFitScouts(ParentPopulation *parentPopulation) {
    static double relativeMinFitness = 0.0;//AlgorithmInfo::instance()->getMinFitnessOfExistingForkingPopulationsRelative() * obf;

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

    vector<ScoutPopulation*> scoutsToRemove;

    // discard child populations whose fitness are less than minfitness
    for ( std::vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        // check if child population's fitness is below minimum fitness, if so discard
        if ((*it).getFitness() < relativeMinFitness) {
            // will remove this population
            scoutsToRemove.push_back( &(*it) );
        }
    }

    // remove populations whose fitness is below min fitness
    for( vector<ScoutPopulation*>::iterator it = scoutsToRemove.begin(); it != scoutsToRemove.end(); it++) {
        parentPopulation->removeScoutPopulation( *it );
    }
}

void Algorithm::GenerationalGA(Population *population, size_t suggestedPopulationSize, Algorithm::CX_TYPE cxType) {
    vector<Individual> newIndividuals;
    int momIndex, dadIndex;

    size_t tournamentSize = Algorithm::TournamentSize;

    Individual bestBefore = population->best();
    population->setPreviousBest( bestBefore );

    // population size may be odd (i+1) is for that reason
    for (int i = 0; i < suggestedPopulationSize && (i+1) < suggestedPopulationSize; i+=2) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            switch (cxType) {
                case HillClimbing:
                    offspring = HillClimbingCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );
                    break;
                default:
                    offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );
            }


            offspring[0].gaussianMutate(population->getMutationStepSize());
            offspring[1].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
            offspring[1].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]) || !population->isValidIndividual(offspring[1]));


        newIndividuals.push_back( offspring[0] );
        newIndividuals.push_back( offspring[1] );
    }

    // population size may be odd, if so init last individual
    if ( suggestedPopulationSize % 2 ) {
        // select parents
        vector<Individual> offspring;
        do {
            SelectParents( population, &momIndex, &dadIndex, tournamentSize );

            switch (CX_Type) {
                case HillClimbing:
                    offspring = HillClimbingCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );
                    break;
                default:
                    offspring = SimulatedBinaryCrossover( population->getIndividual(momIndex), population->getIndividual(dadIndex) );
            }

            offspring[0].gaussianMutate(population->getMutationStepSize());

            offspring[0].setPopulation( population );
        } while(!population->isValidIndividual(offspring[0]));

        newIndividuals.push_back( offspring[0] );
    }

    population->setIndividuals( newIndividuals );

    ApplyElitism( population, bestBefore );

//    assert( newIndividuals.size() == Algorithm::PopulationSize );
}

void Algorithm::SOSAdjustSearchSpace(ParentPopulation *parentPopulation) {
    static double globalMinRadius = Algorithm::SOS_MinDiameterRelative
                                    * ( Algorithm::MaxCoordinate-Algorithm::MinCoordinate );

    static size_t minScoutSize = (size_t) (Algorithm::SOS_MinScoutPopulationSizeRelative * Algorithm::PopulationSize);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for(vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        vector<Individual> individualsToParent;

        // update center
        (*it).setCenter( (*it).best() );

        // get invalid individuals
        vector<Individual>& scoutIndividuals = (*it).getIndividuals();
        for( vector<Individual>::iterator ind = scoutIndividuals.begin(); ind != scoutIndividuals.end(); ind++ ) {
            if ( !(*it).isValidIndividual( *ind ) ) {
                individualsToParent.push_back( *ind );
            }
        }

        assert( (*it).getPopulationSize() >= 4 );
        size_t calculatedScoutSize = (*it).getPopulationSize() - individualsToParent.size();

        // add random individuals if number of individuals in scout dropped below min scout population threshold.
        if ( calculatedScoutSize < minScoutSize ) {
            size_t diff = minScoutSize - calculatedScoutSize;
            for (int i = 0; i < diff; i++) {
                (*it).initIndividualRandom( &(individualsToParent[individualsToParent.size()-i-1]) );
            }
            for (int i = 0; i < diff; i++) {
                individualsToParent.pop_back();
            }
        }

        // remove individuals from scout
        (*it).removeIndividuals( individualsToParent );

        // add these individuals to parent
        parentPopulation->addIndividuals( individualsToParent );
    }
}

double Algorithm::SOSGetRelativeDynamism(Population *population, double totalDynamism) {
    if ( totalDynamism <= 0 )
        return (double) 1./Algorithm::DimensionSize;

    return population->getDynamism()/totalDynamism;
}

double Algorithm::SOSGetRelativeFitness(Population *population, double overallMinFitness, double totalFitness) {
    if ( totalFitness <= 0 )
        return (double) 1./Algorithm::DimensionSize;

    return (population->getFitness()-overallMinFitness)/totalFitness;
}

double Algorithm::SOSGetRelativeQuality(Population *population,
                                        double beta, double totalDynamism, double overallMinFitness,
                                        double totalFitness) {
    return Algorithm::SOS_ALpha * beta * SOSGetRelativeDynamism(population, totalDynamism)
           + (1.-Algorithm::SOS_ALpha*beta)* SOSGetRelativeFitness(population, overallMinFitness, totalFitness);
}

double Algorithm::SOSGetQuality(Population *population, double beta, double totalDynamism, double overallMinFitness,
                                double totalFitness) {
    return Algorithm::SOS_ALpha * beta *  population->getDynamism()
           + (1.-Algorithm::SOS_ALpha*beta)* SOSGetRelativeFitness(population, overallMinFitness, totalFitness);
}

size_t Algorithm::SOSSuggestedSize(Population *population, double quality) {
    return (size_t)((quality*Algorithm::PopulationSize+population->getPopulationSize()) / 2);
}

ScoutPopulation Algorithm::SOSFork(ParentPopulation *parentPopulation) {
    ScoutPopulation scoutPopulation;

    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();

    // get most dense scout population

    // get elements whose fitness is better than min scout population fitness
    double overallBestIndividualFitness = parentPopulation->overallBestFitness();
    double minFitnessOfScout = Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative * overallBestIndividualFitness;
    if (overallBestIndividualFitness < 0) {
        minFitnessOfScout = (1-Algorithm::SOS_MinFitnessOfNewForkingPopulationsRelative) * overallBestIndividualFitness + overallBestIndividualFitness;
    }

    struct FitIndividualInfo {
        int index; // index in m_individuals array
        double radius; // maximum radius calculated
        double maxRadius; // possible maximum radius (this value is used for next possible scout center)
        vector<Individual*> individuals;
    };

    // holds the index of fit individuals
    vector<FitIndividualInfo> fitIndividuals;
    vector<Individual>& individuals = parentPopulation->getIndividuals();
    for ( int i = 0; i < parentPopulation->getPopulationSize(); i++ ) {
        // cout << "curr fitness: " << m_individuals[i].getFitness() << " minfitness: " << minFitnessOfScout << endl;
        if ( individuals[i].getFitness() >= minFitnessOfScout ) {
            FitIndividualInfo fii;
            fii.index = i;
            fii.radius = 0;
            fii.maxRadius = globalMaxDiameter;
            fii.individuals.push_back( &individuals[i] ); // add center to scout
            fitIndividuals.push_back( fii );
        }
    }

    // no suitable individual
    if (fitIndividuals.size() <= 0) {
        return scoutPopulation;
    }

    // fit individuals to other individuals distance array (distance: euclidean)
    double distanceArray[ fitIndividuals.size() ][ individuals.size() ];


    // fill distance array
    // get maximum possible radius
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < fitIndividuals.size(); j++) {
            if ( i == j ) {// distance to itself
                distanceArray[i][j] = 0;
            }
            else {
                distanceArray[i][j] = individuals[fitIndividuals[i].index].getEuclideanDistanceToIndividual(individuals[fitIndividuals[j].index]);

                // compared candidate center individual to other members, if compared member has better
                // fitness radius will not be incremented through this member
                if ( individuals[fitIndividuals[i].index].getFitness() < individuals[fitIndividuals[j].index].getFitness() ) {
                    if ( fitIndividuals[i].maxRadius > distanceArray[i][j] ) {
                        fitIndividuals[i].maxRadius = distanceArray[i][j];
                    }
                }
            }
        }
    }

    // get variables to calculate density
    for (size_t i = 0; i < fitIndividuals.size(); i++) {
        for (size_t j = 0; j < fitIndividuals.size(); j++) {
            if ( i == j ) {// distance to itself
                continue;
            }
            else {
                // FIXME burada ara elemanların eklenmesi density'i olumsuz etkileyebilir. Bu da dikkate alınmalı
                // check max scout size && radius constraints
                if ( distanceArray[i][j] < fitIndividuals[i].maxRadius
                     && fitIndividuals[i].individuals.size() < maxScoutSize) {
                    // update itemcount in scout (used in desity calculation)
                    fitIndividuals[i].individuals.push_back( &(individuals[fitIndividuals[j].index]) );

                    // update radius
                    if ( fitIndividuals[i].radius < distanceArray[i][j] ) {
                        fitIndividuals[i].radius = distanceArray[i][j];
                    }
                }
            }
        }
    }


    // get density
    FitIndividualInfo mostDenseScout;
    double maxDensity = 0;
    bool foundValidScout = false;
    for (vector<FitIndividualInfo>::iterator it = fitIndividuals.begin(); it != fitIndividuals.end(); it++) {
        // check radius constraints
        if ((*it).radius <= 0
            || (*it).radius < globalMinDiameter )
            continue;

        // check number of individuals constraints
        if ((*it).individuals.size() < minScoutSize)
            continue;

        double d = (double) (*it).individuals.size() / (*it).radius;
        if ( d > maxDensity ){
            foundValidScout = true;
            mostDenseScout = *it;
            maxDensity = d;
        }
    }

    // no suitable scout to fork
    if (!foundValidScout) {
        return scoutPopulation;
    }


    scoutPopulation.setIndividuals( mostDenseScout.individuals );
    scoutPopulation.setSuggestedSize( mostDenseScout.individuals.size() );
    scoutPopulation.setDiameter( mostDenseScout.radius );


    vector<Individual> newIndividuals;
    // prepare individuals for parent population (scout now owns the individuals)
    size_t newIndividualsCounter = 0;

    for (int i = 0; i < individuals.size(); i++) {
        bool found = false;
        for (vector<Individual*>::iterator it = mostDenseScout.individuals.begin();
             it != mostDenseScout.individuals.end(); it++) {
            if ( &(individuals[i]) == *it ) {
                found = true;
                break;
            }
        }

        if (found)
            continue;

        newIndividuals.push_back( individuals[i] );
    }

    // update individuals
    parentPopulation->setIndividuals( newIndividuals );

    // check min parent pop size
    while ( parentPopulation->getPopulationSize() < minBaseSize
            || parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() < minBaseSize) {
        // delete worst scout
        parentPopulation->deleteWorstScoutPopulation();
    }

    parentPopulation->setSuggestedSize( parentPopulation->getSuggestedSize()-mostDenseScout.individuals.size() );

    return scoutPopulation;
}

bool Algorithm::SOSMergeScoutPopulations(ParentPopulation *parentPopulation) {
    ScoutPopulation* scout1, * scout2;
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (size_t i = 0; i < scouts.size(); i++) {
        scout1 = &(scouts[i]);
        const Individual& scout1CenterIndividual = scout1->getCenter();
        for (size_t j = i+1; j < scouts.size(); j++) {
            scout2 = &(scouts[j]);
            const Individual& scout2CenterIndividual = scout2->getCenter();

            double scoutToScoutDistance =  scout1CenterIndividual.getEuclideanDistanceToIndividual( scout2CenterIndividual );

            // do we need merge?
            if ( scoutToScoutDistance < scout1->getDiameter() || scoutToScoutDistance < scout2->getDiameter()) {
                if ( scout1->getFitness() < scout2->getFitness() ) { // merge into scout2
                    MergeScoutPopulations( parentPopulation, scout2, scout1 );
                }
                else { // merge into scout1
                    MergeScoutPopulations( parentPopulation, scout1, scout2 );
                }

                return true;
            }
        }
    }

    return false;
}

void Algorithm::MergeScoutPopulations(ParentPopulation* parentPopulation, ScoutPopulation *better, ScoutPopulation *worse) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);
    static double globalMaxDiameter = SOS_MaxDiameterRelative * (MaxCoordinate-MinCoordinate);
    static size_t maxScoutSize = (size_t)(Algorithm::PopulationSize * Algorithm::SOS_MaxScoutPopulationSizeRelative);
    static size_t minScoutSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinScoutPopulationSizeRelative);
    static size_t minBaseSize = (size_t) (Algorithm::PopulationSize * Algorithm::SOS_MinBasePopulationSizeRelative);

    // merge into better scout population

    double newRadius = pow(pow( better->getDiameter(), Algorithm::DimensionSize )
                           + pow( worse->getDiameter(), Algorithm::DimensionSize ), (double)1/Algorithm::DimensionSize);

    newRadius = max(min(newRadius, globalMaxDiameter), globalMinDiameter);

    double distance;
    vector<Individual>& betterScoutIndividuals = better->getIndividuals();
    vector<Individual>& worseScoutIndividuals = worse->getIndividuals();

    better->setDiameter( newRadius );
    better->setSuggestedSize( better->getSuggestedSize() + worse->getSuggestedSize() );

    vector<Individual> individualsToAdd;
    vector<Individual> individualsToParent;
    const Individual& centerInd = better->getCenter();

    // add individuals of better scout
    for (int i = 0; i < better->getPopulationSize(); i++) {
        distance = centerInd.getEuclideanDistanceToIndividual( betterScoutIndividuals[i] );
        if ( distance < better->getDiameter() ) { // add this individual to scout
            individualsToAdd.push_back( betterScoutIndividuals[i] );
        }
        else {
            individualsToParent.push_back( betterScoutIndividuals[i] );
        }
    }

    // add individuals of worse scout
    for (int i = 0; i < worse->getPopulationSize(); i++) {
        distance = centerInd.getEuclideanDistanceToIndividual( worseScoutIndividuals[i] );
        if ( distance < better->getDiameter() ) { // add this individual to scout
            individualsToAdd.push_back( worseScoutIndividuals[i] );
        }
        else {
            individualsToParent.push_back( worseScoutIndividuals[i] );
        }
    }

    better->setIndividuals( individualsToAdd );
    parentPopulation->addIndividuals( individualsToParent );

    // delete worst scout population
    for( vector<ScoutPopulation>::iterator it = parentPopulation->getScoutPopulations().begin();
         it != parentPopulation->getScoutPopulations().end(); it++) {
        if ( &(*it) == worse ) {
            parentPopulation->getScoutPopulations().erase( it );
            break;
        }
    }
}


void Algorithm::SOSReduceScoutDiameters(ParentPopulation *parentPopulation) {
    static double globalMinDiameter = SOS_MinDiameterRelative * (MaxCoordinate-MinCoordinate);

    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); it++) {
        // update radius
        (*it).setDiameter( max((*it).getDiameter()*Algorithm::SOS_DiameterReduceFactor, globalMinDiameter) );
    }
}

std::vector<Individual> Algorithm::SimulatedBinaryCrossover(const Individual &mom, const Individual &dad) {
    vector<Individual> offsprings;
    Individual sis = mom;
    Individual bro = dad;

    if (((double) rand() / RAND_MAX) < Algorithm::CrossoverProbability) {
        size_t lower, upper, rnd1, rnd2;

        double u, beta;

        rnd1 = rand() % Individual::getDimension();
        rnd2 = rand() % Individual::getDimension();

        if (rnd1 >= rnd2) {
            upper = rnd1;
            lower = rnd2;
        } else {
            upper = rnd2;
            lower = rnd1;
        }

        sis.setFitnessCalculated(false);
        bro.setFitnessCalculated(false);

        double* sisValues = sis.getValues();
        double* broValues = bro.getValues();
        double* momValues = mom.getValues();
        double* dadValues = dad.getValues();


        u = (double) rand() / RAND_MAX;
        if (u <= 0.5) {
            beta = pow(2 * u, (double) 1 / 3);
        } else {
            beta = (double) 1 / pow(2 * (1 - u), (double) 1 / 3);
        }
        for (size_t i = lower; i <= upper; i++) {
            sisValues[i] = 0.5 * ((1 + beta) * momValues[i] + (1 - beta) * dadValues[i]);
            broValues[i] = 0.5 * ((1 - beta) * momValues[i] + (1 + beta) * dadValues[i]);

            // fix min max coord values
            if (sisValues[i] < mincoordinate || sisValues[i] > maxcoordinate)
                sisValues[i] = (momValues[i] + dadValues[i]) / 2;

            if (broValues[i] < mincoordinate || broValues[i] > maxcoordinate)
                broValues[i] = (momValues[i] + dadValues[i]) / 2;
        }
    }

    offsprings.push_back( sis );
    offsprings.push_back( bro );

    return offsprings;
}

std::vector<Individual> Algorithm::HillClimbingCrossover(Individual &mom, Individual &dad) {
    int lower, upper, rnd1, rnd2;

    Individual offsprings[10];

    vector<Individual> result;

    Individual oldMom = mom;
    Individual oldDad = dad;

    double u,beta;
    for(int k = 0; k < 2; k++){
        for(int j = 0; j < 10; j = j + 2){
            rnd1 = rand() % Algorithm::DimensionSize;
            rnd2 = rand() % Algorithm::DimensionSize;

            if (rnd1 >= rnd2){
                upper = rnd1;
                lower = rnd2;
            }
            else {
                upper = rnd2;
                lower = rnd1;
            }

            offsprings[j] = oldMom;
            offsprings[j+1] = oldDad;

            double* sisValues = offsprings[j].getValues();
            double* broValues = offsprings[j+1].getValues();
            double* momValues = oldMom.getValues();
            double* dadValues = oldDad.getValues();


            u = (double) rand()/RAND_MAX;
            if (u <= 0.5){
                beta = (double)pow(2*u,(double)1/3);
            }else{
                beta = (double)1/pow(2*(1-u),(double)1/3);
            }
            for(int i = lower; i <= upper; i++){
                sisValues[i] = 0.5 * ((1+beta) * momValues[i] + (1-beta) * dadValues[i]);
                broValues[i] = 0.5 * ((1-beta) * momValues[i] + (1+beta) * dadValues[i]);


                if(sisValues[i] < 0)
                    sisValues[i] = (momValues[i]+dadValues[i])/2;
                if(sisValues[i] > 100)
                    sisValues[i] = (momValues[i]+dadValues[i])/2;
                if(broValues[i] < 0)
                    broValues[i] = (momValues[i]+dadValues[i])/2;
                if(broValues[i] > 100)
                    broValues[i] = (momValues[i]+dadValues[i])/2;
            }

            offsprings[j].updateFitness();
            offsprings[j+1].updateFitness();
        }

        // find best individual
        double bestFitness = -DBL_MAX;
        int bestIndex = -1;
        for (int i = 0; i < 10; i++) {
            if (offsprings[i].getFitness() > bestFitness) {
                bestFitness = offsprings[i].getFitness();
                bestIndex = i;
            }
        }

        assert( bestIndex >= 0 && bestIndex < 10 );

        if( (oldMom.getFitness() < oldDad.getFitness()) && (offsprings[bestIndex].getFitness() > oldMom.getFitness())) {
            oldMom = offsprings[bestIndex];
            oldMom.updateFitness();
        }else if((oldMom.getFitness() >= oldDad.getFitness()) && (offsprings[bestIndex].getFitness() > oldDad.getFitness())) {
            oldDad = offsprings[bestIndex];
            oldDad.updateFitness();
        }
    }

    result.push_back( oldMom );
    result.push_back( oldDad );

    return result;
}


void Algorithm::SOSAddRandomImmigrants(Population *population, size_t size ) {
    for (int i = 0; i < size; ++i) {
        Individual ind;
        ind.setPopulation( population );
        population->initIndividualRandom( &ind );

        population->addIndividual( ind );
    }
}

void Algorithm::SOSRemoveWorstIndividuals(Population *population, size_t size) {
    vector<Individual*> worst = Util::getWorstIndividuals( size, population->getIndividuals() );

    vector<Individual> newIndividuals;

    vector<Individual>& individuals = population->getIndividuals();
    for(int i = 0; i < population->getPopulationSize(); i++) {
        bool found = false;
        for (int j = 0; j < size; j++) {
            if ( &(individuals[i]) == worst[j] ) {
                found = true;
                break;
            }
        }

        if (found)
            continue;

        newIndividuals.push_back( individuals[i] );
    }

    population->removeIndividuals( worst );
}


void Algorithm::RunSelfOrganizingScoutsRandomImmigrants(bool fix) {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    ParentPopulation parentPopulation( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        SOSGenerationalGARandomImmigrants( &parentPopulation );

        if (fix)
            parentPopulation.fixIndividuals();

        SOSReduceScoutDiameters( &parentPopulation );

        SOSAdjustSearchSpace( &parentPopulation );

        if (fix)
            parentPopulation.fixIndividuals();

        SOSDeleteNonFitScouts( &parentPopulation );

        if (fix)
            parentPopulation.fixIndividuals();

        SOSAdjustPopulationSizesRandomImmigrants( &parentPopulation );

        if (fix)
            parentPopulation.fixIndividuals();

        ScoutPopulation sp = SOSFork( &parentPopulation );

        if (sp.getIndividuals().size() > 0) { // forked a child
            parentPopulation.addScoutPopulation( sp );
            parentPopulation.generatePopulationRandom();
        }

        if (fix)
            parentPopulation.fixIndividuals();

        while (SOSMergeScoutPopulations( &parentPopulation )) {
        }

        if (fix)
            parentPopulation.fixIndividuals();

        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) { // env. change
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            // merge explicit mem and memory pop select best n
            change_peaks();

            parentPopulation.updateAllFitnesses();

            SOSAdjustSearchSpace( &parentPopulation );

            if (fix)
                parentPopulation.fixIndividuals();

            statistics.addBestErrorAtChangeStat( parentPopulation.overallBestFitness() );
        }

        assert( parentPopulation.getPopulationSize() >= 8 );

        statistics.addStat( parentPopulation.overallBestFitness() );
        statistics.addAbsoluteRecoveryRateStat( parentPopulation.overallBestFitness() );
//        vector<ScoutPopulation>& scouts = parentPopulation.getScoutPopulations();
//        for(int i = 0; i < scouts.size(); i++) {
//                assert( scouts[i].getPopulationSize() >= 4 &&  scouts[i].getPopulationSize() <= 10);
//        }
    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


//void Algorithm::AssignRanks(Population *population) {
//    int i,rank_index,rank;
//    double temp = -1000000,sum=0;
//
//    vector<Individual>& individuals = population->getIndividuals();
//
//    for(rank = 0;rank<population->getPopulationSize();rank++){
//        for(i=0;i<population->getPopulationSize();i++) {
//            if(individuals[i].getFitness() > temp) {
//                rank_index = i;
//                temp = individuals[i].getFitness();
//            }
//        }
//
//        individuals[rank_index].setRank(rank);
//        individuals[rank_index].setLinRank((((double)2/population->getPopulationSize())-((double)(rank*2)/(population->getPopulationSize()*(population->getPopulationSize()-1)))));
//        temp = -100000;
//    }
//}

void Algorithm::RunHyperMutation() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    // generate population
    Population population( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) {
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            change_peaks();
            population.updateFitnesses();
            Individual::setHyperMutationEnabled( true );

            statistics.addBestErrorAtChangeStat( population.best().getFitness() );
        }

        GenerationalGA( &population );

        Individual::setHyperMutationEnabled( false );

        statistics.addStat( population.best().getFitness() );
        statistics.addAbsoluteRecoveryRateStat( population.best().getFitness() );
    }

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}


void Algorithm::RunSelfOrganizingScoutsRandomImmigrantsMovingCenter() {
    srand((unsigned int)time(0));

    Algorithm::statistics.clear();

    Algorithm::InitMovPeaks();

    ParentPopulation parentPopulation( Algorithm::PopulationSize );

    clock_t start = clock();

    // run generations
    for (int i = 0; i < Algorithm::GenerationCount; i++) {
        SOSGenerationalGARandomImmigrantsMovingCenter( &parentPopulation );

        parentPopulation.fixIndividuals();

        SOSReduceScoutDiameters( &parentPopulation );

        SOSAdjustSearchSpace( &parentPopulation );

        parentPopulation.fixIndividuals();

        SOSDeleteNonFitScouts( &parentPopulation );

        parentPopulation.fixIndividuals();

        SOSAdjustPopulationSizesRandomImmigrants( &parentPopulation );

        parentPopulation.fixIndividuals();

        ScoutPopulation sp = SOSFork( &parentPopulation );

        if (sp.getIndividuals().size() > 0) { // forked a child
            parentPopulation.addScoutPopulation( sp );
            parentPopulation.generatePopulationRandom();
        }

        parentPopulation.fixIndividuals();

        while (SOSMergeScoutPopulations( &parentPopulation )) {
        }

        parentPopulation.fixIndividuals();

        // any change?
        if (i > 0 && (i % Algorithm::ChangePeriod == 0)) { // env. change
            statistics.addCalculatedAbsoluteRecoveryRate( statistics.getAbsoluteRecoveryRate() );
            statistics.clearAbsoluteRecoveryRateStat();

            // merge explicit mem and memory pop select best n
            change_peaks();

            parentPopulation.updateAllFitnesses();

            SOSAdjustSearchSpace( &parentPopulation );

            parentPopulation.fixIndividuals();

            Algorithm::statistics.addBestErrorAtChangeStat( parentPopulation.overallBestFitness() );
        }

        assert( parentPopulation.getPopulationSize() >= 8 );

        statistics.addStat( parentPopulation.overallBestFitness() );
        statistics.addAbsoluteRecoveryRateStat( parentPopulation.overallBestFitness() );

    };

    statistics.addTimeSpan( start );

    statistics.addRunOfflineError( statistics.getOfflineError() );
    statistics.addRunAccuracy( statistics.getAccuracy() );
    statistics.addRunAbsoluteRecoveryRate( statistics.getCalculatedAbsoluteRecoveryRate() );
    statistics.addRunBestErrorAtChange( statistics.getBestErrorAtChange() );

    Algorithm::ReleaseMovPeaks();
}

void Algorithm::SOSGenerationalGARandomImmigrantsMovingCenter(ParentPopulation *parentPopulation) {
// generational ga steps for
    vector<ScoutPopulation>& scouts = parentPopulation->getScoutPopulations();
    for (vector<ScoutPopulation>::iterator it = scouts.begin(); it != scouts.end(); ++it) {
        GenerationalGAMovingCenter( &(*it)  );
    }

    GenerationalGAMovingCenter( parentPopulation );
}
