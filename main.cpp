// cd C:\Users\enoch\OneDrive\Documents\ScramblerCalibrator
// .\ScramblerCalibrator.exe
#include "ScramblerCalibrator.h"
#include <iostream>
#include <limits>

int main() {
    ScramblerCalibrator calibrator;
    int choice = 0;

    do {
        std::cout << "\n=== Scrambler Menu ===\n";
        std::cout << "1) Add run\n";
        std::cout << "2) View runs\n";
        std::cout << "3) Compute recommendation\n";
        std::cout << "4) Clear all past data\n";
        std::cout << "5) Exit\n";
        std::cout << "Choice: ";

        std::cin >> choice;

        // Validate input
        if (std::cin.fail()) {
            std::cin.clear(); // clear error flag
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // discard bad input
            std::cout << "Invalid input. Please enter a number between 1 and 5.\n";
            continue;
        }

        switch(choice) {
            case 1: {
                double rotations, errorCm;
                std::cout << "Enter rotations: ";
                std::cin >> rotations;
                std::cout << "Enter error (cm): ";
                std::cin >> errorCm;
                calibrator.addRun(rotations, errorCm);
                break;
            }
            case 2:
                calibrator.viewRuns();
                break;
            case 3:
                calibrator.computeRecommendation();
                break;
            case 4:
                calibrator.clearAllData();
                break;
            case 5:
                std::cout << "Exiting program...\n";
                break;
            default:
                std::cout << "Invalid choice. Please enter a number between 1 and 5.\n";
        }

    } while(choice != 5);

    return 0;
}