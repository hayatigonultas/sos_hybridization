#ifndef EVOLUTIONARYLIB_STATISTICS_H
#define EVOLUTIONARYLIB_STATISTICS_H

#include <stddef.h>
#include <ctime>

#include <vector>

class Statistics {
public:
    Statistics();
    ~Statistics();

    double getAbsoluteRecoveryRate() const;
    double getRunAbsoluteRecoveryRate() const;
    double getCalculatedAbsoluteRecoveryRate() const;

    double getAccuracy() const;
    double getRunAccuracy() const;

    double getBestErrorAtChange() const;
    double getRunBestErrorAtChange() const;

    double getOfflineError() const;
    double getRunOfflineError() const;

    double getOfflinePerformance() const;
    double getRunOfflinePerformance() const;

    double getRunTimespan() const;

    void addStat(double bestFitness);
    void addAbsoluteRecoveryRateStat(double bestFitness);
    void addBestErrorAtChangeStat( double bestFitness );
    void addCalculatedAbsoluteRecoveryRate(double arr);

    void clear();
    void clearAbsoluteRecoveryRateStat();
    void clearRunStat();

    void addRunOfflinePerformance(double performance);
    void addRunOfflineError(double error);
    void addRunAccuracy(double error);
    void addRunBestErrorAtChange(double error);
    void addRunAbsoluteRecoveryRate( double val );

    void addTimeSpan( clock_t start );

private:
    std::vector<double> bestErrorAtChange;
    std::vector<double> accuracy;
    std::vector<double> offlineError;
    std::vector<double> offlinePerformance;
    std::vector<double> absoluteRecoveryRate;
    std::vector<double> calculatedAbsoluteRecoveryRate;

    std::vector<double> runBestErrorAtChange;
    std::vector<double> runOfflineError;
    std::vector<double> runOfflinePerformance;
    std::vector<double> runAccuracy;
    std::vector<double> runAbsoluteRecoveryRate;

    std::vector<double> runTimeSpan;


    size_t m_generationCounter;
};


#endif //EVOLUTIONARYLIB_STATISTICS_H
