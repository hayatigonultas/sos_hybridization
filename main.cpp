#include <iostream>

#include "Algorithm.h"
#include "Util.h"
#include "Test.h"

using namespace std;

static const char* filename = "mpbrc";

#ifndef TEST

int main() {
    Util::Load_RC_File( filename );

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSimpleGA();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SGA OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;


    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunHyperMutation();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "HM OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunRandomImmigrants();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "RI OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunMemorySearch();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "MS OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScouts();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SOS OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

    Algorithm::CX_Type = Algorithm::HillClimbing;
    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScouts();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SOS+LS OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;


    Algorithm::CX_Type = Algorithm::SBX;
    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScoutsRandomImmigrants(false);
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SOS+RI+NO_FIX OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;


    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScoutsRandomImmigrants(true);
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SOS+RI+FIX OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;

    Algorithm::statistics.clearRunStat();
    for (size_t i = 0; i < Algorithm::RunCount; i++) {
        Algorithm::RunSelfOrganizingScoutsRandomImmigrantsMovingCenter();
        Util::LoadBar( i, Algorithm::RunCount );
    }

    cout << "SOS+RI+MOVING_CENTER_FIX OE: " << Algorithm::statistics.getRunOfflineError()
    << " AC: " << Algorithm::statistics.getRunAccuracy()
    << " BEAC: " << Algorithm::statistics.getRunBestErrorAtChange()
    << " ARR: " << Algorithm::statistics.getRunAbsoluteRecoveryRate()
    << " RT: " << Algorithm::statistics.getRunTimespan() <<  endl;


    return 0;
}

#endif
