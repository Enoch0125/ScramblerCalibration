#ifndef SCRAMBLERCALIBRATOR_H
#define SCRAMBLERCALIBRATOR_H

#include <vector>
#include <string>

struct RunRecord {
    double rotations;
    double errorCm;
    RunRecord(double r = 0.0, double e = 0.0) : rotations(r), errorCm(e) {}
};

class ScramblerCalibrator {
public:
    ScramblerCalibrator();

    void loadPreviousData();  // load run_data.csv if exists
    void addRun(double rotations, double errorCm);
    void viewRuns() const;
    void computeRecommendation(); // IRLS + delta + logs
    void clearAllData();
private:
    std::vector<RunRecord> runs;

    // --- Constants ---
    const double WHEEL_DIAMETER_IN = 2.875;
    const double IN_TO_CM = 2.54;
    const double WHEEL_CIRCUMFERENCE_CM;
    const double TUKEY_C = 4.685;
    const double CONVERGENCE_TOL = 1e-6;
    const int MAX_ITER = 100;
    const double CHANGE_LIMIT = 2.0; // +/- rotations warning

    // --- Helper functions ---
    std::pair<double,double> ordinaryLeastSquares() const;
    std::pair<double,double> weightedLeastSquares(const std::vector<double>& weights) const;
    double computeMAD(const std::vector<double>& residuals) const;
    double fallbackStdEstimate(const std::vector<double>& residuals) const;
    void appendRunDataCsv(const RunRecord& r) const;
    void appendIterationLog(int regressionId, int iteration, double slope, double intercept,
                            double avgAbsResidual, int weightedPoints) const;
    int estimateNextRegressionId() const;
};

#endif // SCRAMBLERCALIBRATOR_H