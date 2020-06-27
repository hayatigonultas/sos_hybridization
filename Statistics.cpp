#include "Statistics.h"

extern double global_max;

using namespace std;

Statistics::Statistics() {
}

Statistics::~Statistics() {
}


double Statistics::getOfflineError() const {
    double d = 0.;
    for (vector<double>::const_iterator it = offlineError.begin();
         it != offlineError.end(); it++) {
        d += *it;
    }

    return d / (offlineError.size());
}

double Statistics::getOfflinePerformance() const {
    double d = 0.;

    for (vector<double>::const_iterator it = offlinePerformance.begin();
         it != offlinePerformance.end(); it++) {
        d += *it;
    }

    return d / (offlinePerformance.size());
}


void Statistics::addStat(double bestFitness) {
    offlineError.push_back( global_max - bestFitness );
    offlinePerformance.push_back( bestFitness );
    accuracy.push_back( bestFitness/global_max );
}

void Statistics::clear() {
    offlineError.clear();
    offlinePerformance.clear();
    accuracy.clear();
    bestErrorAtChange.clear();
    calculatedAbsoluteRecoveryRate.clear();
    absoluteRecoveryRate.clear();


    runTimeSpan.clear();
}

void Statistics::addRunOfflinePerformance(double performance) {
    runOfflinePerformance.push_back( performance );
}


void Statistics::addRunOfflineError(double error) {
    runOfflineError.push_back( error );
}

double Statistics::getRunOfflineError() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runOfflineError.begin();
         it != runOfflineError.end(); it++) {
        d += *it;
    }

    return d / (runOfflineError.size());
}


double Statistics::getRunOfflinePerformance() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runOfflinePerformance.begin();
         it != runOfflinePerformance.end(); it++) {
        d += *it;
    }

    return d / (runOfflinePerformance.size());
}

void Statistics::clearRunStat() {
    runOfflinePerformance.clear();
    runOfflineError.clear();
    runBestErrorAtChange.clear();
    runAccuracy.clear();
    runAbsoluteRecoveryRate.clear();

    runTimeSpan.clear();
}

double Statistics::getAccuracy() const {
    double a = 0.;
    for (vector<double>::const_iterator it = accuracy.begin(); it != accuracy.end(); it++) {
        a += *it;
    }

    return a / (accuracy.size());
}

double Statistics::getRunAccuracy() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runAccuracy.begin(); it != runAccuracy.end(); it++) {
        a += *it;
    }

    return a / (runAccuracy.size());
}

void Statistics::addRunAccuracy(double error) {
    runAccuracy.push_back( error );
}

void Statistics::addTimeSpan(clock_t start) {
    clock_t end = clock();
    double diffTime = 0.;

    diffTime = ((double)(end - start))/CLOCKS_PER_SEC;
    runTimeSpan.push_back( diffTime );
}

double Statistics::getRunTimespan() const {
    double d = 0.;

    for (vector<double>::const_iterator it = runTimeSpan.begin(); it != runTimeSpan.end(); it++) {
        d += *it;
    }

    return d / (runTimeSpan.size());
}

void Statistics::addRunBestErrorAtChange(double error) {
    runBestErrorAtChange.push_back( error );
}

double Statistics::getBestErrorAtChange() const {
    double a = 0.;
    for (vector<double>::const_iterator it = bestErrorAtChange.begin(); it != bestErrorAtChange.end(); it++) {
        a += *it;
    }

    return a / (bestErrorAtChange.size());
}

double Statistics::getRunBestErrorAtChange() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runBestErrorAtChange.begin(); it != runBestErrorAtChange.end(); it++) {
        a += *it;
    }

    return a / (runBestErrorAtChange.size());
}

void Statistics::addBestErrorAtChangeStat(double bestFitness) {
    bestErrorAtChange.push_back( global_max - bestFitness );
}

double Statistics::getAbsoluteRecoveryRate() const {
    double bestAtGenZero = absoluteRecoveryRate.at(0);

    double a = 0.;
    for (vector<double>::const_iterator it = absoluteRecoveryRate.begin(); it != absoluteRecoveryRate.end(); it++) {
        a += *it - bestAtGenZero;
    }

    double b = absoluteRecoveryRate.size() * ( global_max - bestAtGenZero );

    return a / b;
}



double Statistics::getRunAbsoluteRecoveryRate() const {
    double a = 0.;

    for (vector<double>::const_iterator it = runAbsoluteRecoveryRate.begin(); it != runAbsoluteRecoveryRate.end(); it++) {
        a += *it;
    }

    return a / (runAbsoluteRecoveryRate.size());
}

void Statistics::addRunAbsoluteRecoveryRate(double val) {
    runAbsoluteRecoveryRate.push_back( val );
}

void Statistics::addAbsoluteRecoveryRateStat(double bestFitness) {
    absoluteRecoveryRate.push_back( bestFitness );
}

void Statistics::addCalculatedAbsoluteRecoveryRate(double arr) {
    calculatedAbsoluteRecoveryRate.push_back(arr);
}

double Statistics::getCalculatedAbsoluteRecoveryRate() const {
    double a = 0.;

    for (vector<double>::const_iterator it = calculatedAbsoluteRecoveryRate.begin(); it != calculatedAbsoluteRecoveryRate.end(); it++) {
        a += *it;
    }

    return a / (calculatedAbsoluteRecoveryRate.size());
}

void Statistics::clearAbsoluteRecoveryRateStat() {
    absoluteRecoveryRate.clear();
}
