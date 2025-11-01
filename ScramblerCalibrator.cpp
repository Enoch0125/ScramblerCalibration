#include "ScramblerCalibrator.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ScramblerCalibrator::ScramblerCalibrator()
    : WHEEL_CIRCUMFERENCE_CM(M_PI * WHEEL_DIAMETER_IN * IN_TO_CM)
{
    loadPreviousData();
}

// ---------- Load previous runs ----------
void ScramblerCalibrator::loadPreviousData() {
    std::ifstream in("run_data.csv");
    if (!in.is_open()) return;
    std::string line;
    std::getline(in,line); // skip header
    double r,e,d;
    std::string ts;
    while (in >> r) {
        if (!(in >> std::ws)) break; // safety
        in.ignore(1, ',');
        if (!(in >> e)) break;
        in.ignore(1, ',');
        if (!(in >> d)) break;
        in.ignore(1, ',');
        if (!std::getline(in, ts)) break;
        runs.push_back(RunRecord(r,e));
    }
    std::cout << "Loaded " << runs.size() << " previous runs.\n";
}

// ---------- Add new run ----------
void ScramblerCalibrator::addRun(double rotations, double errorCm) {
    RunRecord r(rotations,errorCm);
    runs.push_back(r);
    appendRunDataCsv(r);
    std::cout << "Recorded: rot=" << std::fixed << std::setprecision(2)
              << rotations << ", err=" << errorCm << " cm\n";
}

// ---------- View runs ----------
void ScramblerCalibrator::viewRuns() const {
    if (runs.empty()) { std::cout << "No runs recorded.\n"; return; }
    for (size_t i=0;i<runs.size();i++)
        std::cout << i+1 << ") rot=" << std::fixed << std::setprecision(2)
                  << runs[i].rotations << ", err=" << runs[i].errorCm << " cm\n";
}

// ---------- IRLS Recommendation ----------
void ScramblerCalibrator::computeRecommendation() {
    if (runs.size() < 2) {
        std::cout << "Need at least 2 runs.\n"; return;
    }
    int regressionId = estimateNextRegressionId();

    // Start OLS
    auto [slope, intercept] = ordinaryLeastSquares();
    bool converged = false;

    for (int iter=1; iter<=MAX_ITER; iter++) {
        std::vector<double> residuals;
        for (auto& r: runs) residuals.push_back(r.errorCm - (slope*r.rotations + intercept));

        double mad = computeMAD(residuals);
        if (mad < 1e-9) mad = fallbackStdEstimate(residuals);
        if (mad < 1e-9) mad = 1.0;

        // Tukey weights
        std::vector<double> weights(residuals.size(),0.0);
        double cMad = TUKEY_C*mad;
        int weightedPoints = 0;
        double avgAbsResidual = 0;
        for (size_t i=0;i<residuals.size();i++){
            double u = residuals[i]/cMad;
            if (std::abs(u)<1.0) {
                double tmp = 1-u*u;
                weights[i] = tmp*tmp;
                weightedPoints++;
            } else weights[i]=0.0;
            avgAbsResidual += std::abs(residuals[i]);
        }
        avgAbsResidual /= residuals.size();

        // Weighted least squares
        auto [newSlope,newIntercept] = weightedLeastSquares(weights);

        appendIterationLog(regressionId, iter, newSlope,newIntercept,avgAbsResidual,weightedPoints);

        double change = std::abs(newSlope-slope)+std::abs(newIntercept-intercept);
        slope=newSlope; intercept=newIntercept;
        if (change < CONVERGENCE_TOL) { converged=true; break; }
    }

    if (!converged) std::cout << "Warning: IRLS max iterations reached.\n";

    std::cout << "Final slope=" << slope << ", intercept=" << intercept << "\n";

    if (std::abs(slope)<1e-12) {
        std::cout << "Slope near zero; cannot compute recommendation.\n"; return;
    }

    double rotationsZero = -intercept/slope;
    double delta = rotationsZero - runs.back().rotations;

    std::cout << "Recommended total: " << std::fixed << std::setprecision(2) << rotationsZero << " rotations\n";
    std::cout << "Change: " << std::showpos << std::fixed << std::setprecision(2) << delta
              << " rotations (from last run " << runs.back().rotations << ")\n";
    if (std::abs(delta) > CHANGE_LIMIT)
        std::cout << "WARNING: change exceeds ±" << CHANGE_LIMIT << " rotations!\n";

    double recDistCm = rotationsZero*WHEEL_CIRCUMFERENCE_CM;
    std::cout << "Recommended travel distance: " << recDistCm << " cm\n";
}

// ---------- Math helpers ----------
std::pair<double,double> ScramblerCalibrator::ordinaryLeastSquares() const {
    double sumX=0,sumY=0,sumXY=0,sumX2=0;
    for (auto& r: runs){
        sumX+=r.rotations; sumY+=r.errorCm;
        sumXY+=r.rotations*r.errorCm;
        sumX2+=r.rotations*r.rotations;
    }
    double n = runs.size();
    double denom=n*sumX2-sumX*sumX;
    if (std::abs(denom)<1e-12) return {0.0,sumY/n};
    double slope=(n*sumXY-sumX*sumY)/denom;
    double intercept=(sumY-slope*sumX)/n;
    return {slope,intercept};
}

std::pair<double,double> ScramblerCalibrator::weightedLeastSquares(const std::vector<double>& weights) const {
    double sumW=0,sumWX=0,sumWY=0;
    for (size_t i=0;i<runs.size();i++){
        sumW+=weights[i]; sumWX+=weights[i]*runs[i].rotations; sumWY+=weights[i]*runs[i].errorCm;
    }
    if (sumW<1e-12) return ordinaryLeastSquares();
    double xwMean=sumWX/sumW, ywMean=sumWY/sumW;
    double num=0,den=0;
    for (size_t i=0;i<runs.size();i++){
        num+=weights[i]*(runs[i].rotations-xwMean)*(runs[i].errorCm-ywMean);
        den+=weights[i]*(runs[i].rotations-xwMean)*(runs[i].rotations-xwMean);
    }
    if (std::abs(den)<1e-12) return {0.0,ywMean};
    double slope=num/den;
    double intercept=ywMean-slope*xwMean;
    return {slope,intercept};
}

double ScramblerCalibrator::computeMAD(const std::vector<double>& residuals) const {
    std::vector<double> absRes=residuals;
    for(auto &v:absRes) v=std::abs(v);
    size_t n=absRes.size();
    std::nth_element(absRes.begin(),absRes.begin()+n/2,absRes.end());
    double median=absRes[n/2];
    if(n%2==0){
        std::nth_element(absRes.begin(),absRes.begin()+n/2-1,absRes.end());
        median=(median+absRes[n/2-1])/2.0;
    }
    return median;
}

double ScramblerCalibrator::fallbackStdEstimate(const std::vector<double>& residuals) const {
    double sum=0;
    for(auto r:residuals) sum+=std::abs(r);
    return sum/residuals.size();
}

// ---------- CSV Logging ----------
void ScramblerCalibrator::appendRunDataCsv(const RunRecord& r) const {
    std::ofstream out("run_data.csv", std::ios::app);
    if(out.tellp()==0) out << "rotations,error_cm,distance_cm,timestamp\n"; // header
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    double distance = r.rotations*WHEEL_CIRCUMFERENCE_CM;
    out << r.rotations << "," << r.errorCm << "," << distance << ","
        << std::ctime(&now);
    out.close();
}

void ScramblerCalibrator::appendIterationLog(int regressionId, int iteration, double slope, double intercept,
                                              double avgAbsResidual, int weightedPoints) const {
    std::ofstream out("iteration_log.csv", std::ios::app);
    if(out.tellp()==0) out << "timestamp,regression_id,iteration,slope,intercept,avg_abs_residual,weighted_points\n";
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    out << std::ctime(&now) << "," << regressionId << "," << iteration << ","
        << slope << "," << intercept << "," << avgAbsResidual << "," << weightedPoints << "\n";
    out.close();
}

// ---------- Regression ID ----------
int ScramblerCalibrator::estimateNextRegressionId() const {
    std::ifstream in("iteration_log.csv");
    int maxId=0;
    if(!in.is_open()) return 1;
    std::string line;
    std::getline(in,line); // skip header
    while(std::getline(in,line)){
        std::stringstream ss(line);
        std::string ts;
        int rid;
        if(std::getline(ss, ts, ',')){
            if(ss >> rid) if(rid>maxId) maxId=rid;
        }
    }
    return maxId+1;
    
}

void ScramblerCalibrator::clearAllData() {
    // Clear CSV files
    std::ofstream("run_data.csv", std::ios::trunc);
    std::ofstream("iteration_log.csv", std::ios::trunc);

    // Clear in-memory runs
    runs.clear();

    std::cout << "All past data cleared.\n";
}