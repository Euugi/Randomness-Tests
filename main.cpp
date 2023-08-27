#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

static void print_matrix(const std::vector<std::vector<double>>& t) {
    for (const auto& row : t) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

class BinaryMatrix {

public:
    BinaryMatrix(const std::vector<std::vector<double>>& matrix, int rows, int columns)
            : numberOfRows(rows), numberOfColumns(columns), M(matrix), m(std::min(rows, columns)) {}

    int computeRank() {
        for (int i = 0; i < this->m - 1; i++) {
            if (this->M[i][i] == 1) {
                this->applyRowOperations(i, true);
            } else {
                if (this->findRowToSwapWithUnitElement(i, true) == 1) {
                    this->applyRowOperations(i, true);
                }
            }
        }

        for (int i = this->m - 1; i > 0; i--) {
            if (this->M[i][i] == 1) {
                this->applyRowOperations(i, false);
            } else {
                if (this->findRowToSwapWithUnitElement(i, false) == 1) {
                    this->applyRowOperations(i, false);
                }
            }
        }

        return this->determineRank();
    }

    void applyRowOperations(int i, bool forwardElimination) {
        if (forwardElimination) {
            for (int j = i + 1; j < this->numberOfRows; j++) {
                if (this->M[j][i] == 1) {
                    for (int k = 0; k < this->numberOfColumns; ++k) {
                        this->M[j][k] = int(this->M[j][k] + this->M[i][k]) % 2;
                    }
                }
            }
        } else {
            for (int j = i - 1; j >= 0; j--) {
                if (this->M[j][i] == 1) {
                    for (int k = 0; k < this->numberOfColumns; ++k) {
                        this->M[j][k] = int(this->M[j][k] + this->M[i][k]) % 2;
                    }
                }
            }
        }
    }

    int findRowToSwapWithUnitElement(int rowIndex, bool forwardElimination) {
        int isRowSwapped = 0;
        int index;

        if (forwardElimination) {
            index = rowIndex + 1;
            while (index < this->numberOfRows && this->M[index][rowIndex] == 0) {
                index++;
            }
            if (index < this->numberOfRows) {
                isRowSwapped = this->swapRows(rowIndex, index);
            }
        } else {
            index = rowIndex - 1;
            while (index >= 0 && this->M[index][rowIndex] == 0) {
                index--;
            }
            if (index >= 0) {
                isRowSwapped = this->swapRows(rowIndex, index);
            }
        }
        return isRowSwapped;
    }

    int swapRows(int firstRow, int secondRow) {
        std::vector<double> temp = this->M[firstRow];
        this->M[firstRow] = this->M[secondRow];
        this->M[secondRow] = temp;
        return 1;
    }

    int determineRank() {
        int rank = this->m;
        for(int i = 0; i < this->numberOfRows; i++){
            int all_zeros = 1;
            for (int j = 0; j < this->numberOfColumns; ++j) {
                if (this->M[i][j] == 1) {
                    all_zeros = 0;
                    continue;
                }
            }
            if (all_zeros == 1) {
                rank -= 1;
            }
        }
        return rank;
    }

private:
    int numberOfRows;
    int numberOfColumns;
    std::vector<std::vector<double>> M;
    int m = std::min(numberOfRows, numberOfColumns);
};

int BinaryMatrixRankTest() {
    std::cout << "Binary matrix rank test" << std::endl;


    std::ifstream inputFile("/home/patryk/CLionProjects/untitled9/ciphersBinary.txt");
    if (!inputFile) {
        std::cerr << "Nie można otworzyć pliku wejściowego." << std::endl;
        return 1;
    }

    std::ofstream outputFile("/home/patryk/CLionProjects/BinaryMatrixRankTest/BinaryMatrixRankTest.txt");
    if (!outputFile) {
        std::cerr << "Nie można otworzyć pliku wyjściowego." << std::endl;
        return 1;
    }

    double chi2 = 0;


    std::string line;
    while (std::getline(inputFile, line)) {
        double result = 0.0;
        int matrixDimension = 32;
        std::string binaryData = line;
        std::string::size_type binaryDataLength = binaryData.length();
        int subMatrixSize = matrixDimension * matrixDimension;
        int numberOfSubMatrices = floor(binaryDataLength / subMatrixSize);
        std::string::size_type subMatrixStart = 0, subMatrixEnd = subMatrixSize;

        if (numberOfSubMatrices > 0) {
            std::vector<int> max_ranks(3, 0);
            for (int i = 0; i < numberOfSubMatrices; ++i) {
                std::string subMatrixData = binaryData.substr(subMatrixStart, subMatrixEnd - subMatrixStart);
                std::vector<double> block(subMatrixData.length(), 0.0);
                std::vector<std::vector<double>> m(matrixDimension, std::vector<double>(matrixDimension, 0.0));
                int dataIndex = 0;
                for (int j = 0; j < matrixDimension; ++j) {
                    for (int k = 0; k < matrixDimension; ++k) {
                        dataIndex = j * matrixDimension + k;
                        if (subMatrixData[dataIndex] == '1') {
                            m[j][k] = 1.0;
                        }
                    }
                }

                print_matrix(m);

                BinaryMatrix ranker(m, matrixDimension, matrixDimension);
                int rank = ranker.computeRank();

                if (rank == matrixDimension) {
                    max_ranks[0]++;
                } else if (rank == (matrixDimension - 1)) {
                    max_ranks[1]++;
                } else {
                    max_ranks[2]++;
                }

                subMatrixStart += subMatrixSize;
                subMatrixEnd += subMatrixSize;
            }
            std::vector<double> piks(3);
            piks[0] = 1.0;
            piks[1] = 0.0;
            piks[2] = 0.0;
            for (int x = 1; x < 50; x++) {
                piks[0] *= 1 - (1.0 / static_cast<double>(std::pow(2, x)));
            }
            piks[1] = 2 * piks[0];
            piks[2] = 1 - piks[0] - piks[1];
            double chi = 0.0;
            for (int i = 0; i < piks.size(); ++i) {
                chi += std::pow((max_ranks[i] - piks[i] * numberOfSubMatrices), 2.0) / (piks[i] * numberOfSubMatrices);
                chi2 = chi;
            }
            double p_val = std::exp(-chi / 2);
            result = p_val;
        } else {
            result = -1.0;
        }

        std::string random = "TRUE";
        if(result < 0.01){
            random = "FALSE";
        }


        outputFile<< chi2 << ";&;" << result << " ;&; " << random << " ; " << "\\\\" << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Operacja zakończona sukcesem." << std::endl;



    return 0;
}


int FrequencyMonobitsTest()
{
    std::ifstream inputFile("/home/patryk/CLionProjects/RandomnessTests/ciphersBinary.txt");
    if (!inputFile) {
        std::cerr << "Nie można otworzyć pliku wejściowego." << std::endl;
        return 1;
    }

    std::ofstream outputFile("/home/patryk/CLionProjects/RandomnessTests/frequencyMonobitsTest.txt");
    if (!outputFile) {
        std::cerr << "Nie można otworzyć pliku wyjściowego." << std::endl;
        return 1;
    }
    outputFile << "Frequency monobits test" << std::endl;
    std::string line;
    while (std::getline(inputFile, line)) {
        int numberOfZeros = 0;
        int numberOfOnes = 0;
        int numberOfBits = 0;
        for (char bit : line) {
            if(bit == '0'){
                numberOfZeros++;
                numberOfBits++;
            } else if(bit == '1'){
                numberOfOnes++;
                numberOfBits++;
            }
        }

        int Sn = numberOfOnes - numberOfZeros;
        double sqrt_n = sqrt(numberOfBits);
        double Sobs = fabs(Sn) / sqrt_n;

        double normalDistribution = 0.5 * (1.0 + erf(Sobs / sqrt(2.0)));
        double pValue = 2.0 * (1.0 - normalDistribution);
        std::string random = "false";
        if(pValue>=0.01){
            random = "true";
        }
        outputFile << numberOfZeros << ";" << numberOfOnes << ";" << numberOfBits << ";" << Sn << ";" << Sobs << ";" << pValue << ";" << random << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Operacja zakończona sukcesem." << std::endl;
}

int RunsTest()
{
    std::ifstream inputFile("/home/patryk/CLionProjects/untitled6/randomness-tests/ciphersBinary.txt");
    if (!inputFile) {
        std::cerr << "Nie można otworzyć pliku wejściowego." << std::endl;
        return 1;
    }

    std::ofstream outputFile("/home/patryk/CLionProjects/untitled6/randomness-tests/RunsTest.txt");
    if (!outputFile) {
        std::cerr << "Nie można otworzyć pliku wyjściowego." << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        double numberOfOnes = 0;
        for (char bit : line) {
            if(bit == '1'){
                numberOfOnes++;
            }
        }
        int runs = 1;
        char lastOne = line[0];
        for(int i = 1; i < line.length(); i++){
            if(lastOne != line[i]){
                runs++;
                lastOne = line[i];
            }
        }
        double Pi = numberOfOnes/line.length();

        double erfc = abs(runs-(2.0*line.length()*Pi*(1-Pi)))/(2*(sqrt(2*line.length())*Pi*(1-Pi)));
        double normalDistribution = 0.5 * (1.0 + erf(erfc));
        double pValue = 2.0 * (1.0 - normalDistribution);

        std::string random = "TRUE";
        if(pValue < 0.01){
            random = "FALSE";
        }

        outputFile << numberOfOnes << ";&;" << runs << ";&;" << Pi << ";&;" << pValue << ";&;" << random << ";\\\\" << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Operacja zakończona sukcesem." << std::endl;
}

int FrequencyTestWithinABlock()
{
    std::ifstream inputFile("/home/patryk/CLionProjects/untitled7/ciphersBinary.txt");
    if (!inputFile) {
        std::cerr << "Nie można otworzyć pliku wejściowego." << std::endl;
        return 1;
    }

    std::ofstream outputFile("/home/patryk/CLionProjects/untitled7/FrequencyBlockTestPhotonBeetle.txt");
    if (!outputFile) {
        std::cerr << "Nie można otworzyć pliku wyjściowego." << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        int numberOfCharacter = 0;
        double numberOfOnes = 0;
        double loopPi = 0.0;
        double loopx2obs = 0.0;
        double x2obs = 0.0;

        for (char bit : line) {
            numberOfCharacter++;
            if(bit == '1'){
                numberOfOnes++;
            }

            if(numberOfCharacter == 10000){
                numberOfCharacter = 0;
                loopPi = numberOfOnes/10000;
                loopx2obs = loopx2obs + pow((loopPi-0.5), 2);
                numberOfOnes = 0;
            }
        }
        x2obs = (4 * 10000) * loopx2obs;

        double result = gsl_sf_gamma_inc_Q((40064.0/10000.0)/2.0, x2obs/2);

        std::string random = "TRUE";
        if(result<0.01){
            random = "FALSE";
        }

        outputFile << std::fixed;
        outputFile << (40064.0/10000.0)/2.0 << "; " << x2obs << " ;" << result <<  " ;" + random << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Operacja zakończona sukcesem." << std::endl;
}

int LongestRunOfOnesInABlock() {
    std::ifstream inputFile("/home/patryk/CLionProjects/untitled7/ciphersBinary.txt");
    if (!inputFile) {
        std::cerr << "Nie można otworzyć pliku wejściowego." << std::endl;
        return 1;
    }

    std::ofstream outputFile("/home/patryk/CLionProjects/untitled7/LongestRunOfOnesInABlockTest.txt");
    if (!outputFile) {
        std::cerr << "Nie można otworzyć pliku wyjściowego." << std::endl;
        return 1;
    }

    double n = 40064.0;
    double M = 128.0;
    double Pi0 = 0.1174;
    double Pi1 = 0.2430;
    double Pi2 = 0.2493;
    double Pi3 = 0.1752;
    double Pi4 = 0.1027;
    double Pi5 = 0.1124;
    std::string line;
    while (std::getline(inputFile, line)) {
        double number = 0;
        double v0 = 0;
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        double v4 = 0;
        double v5 = 0;

        double currentRun = 0;
        double maxRun = 0;
        for (int i = 0; i < line.length(); i++) {
            if (line[i] == '1') {
                currentRun++;
                if (currentRun > maxRun) {
                    maxRun = currentRun;
                }
                number++;
            } else {
                currentRun = 0;
                number++;
            }
            if (number == 128) {
                currentRun = 0;
                number = 0;
                if (maxRun <= 4) {
                    v0++;
                }
                if (maxRun == 5) {
                    v1++;
                }
                if (maxRun == 6) {
                    v2++;
                }
                if (maxRun == 7) {
                    v3++;
                }
                if (maxRun == 8) {
                    v4++;
                }
                if (maxRun >= 9) {
                    v5++;
                }
                maxRun = 0;
            }
        }

        double v0obs = pow((v0 - (313 * Pi0)), 2) / (313 * Pi0);
        double v1obs = pow((v1 - (313 * Pi1)), 2) / (313 * Pi1);
        double v2obs = pow((v2 - (313 * Pi2)), 2) / (313 * Pi2);
        double v3obs = pow((v3 - (313 * Pi3)), 2) / (313 * Pi3);
        double v4obs = pow((v4 - (313 * Pi4)), 2) / (313 * Pi4);
        double v5obs = pow((v5 - (313 * Pi5)), 2) / (313 * Pi5);

        double x2obs = v0obs + v1obs + v2obs + v3obs + v4obs + v5obs;

        double pValue = gsl_sf_gamma_inc_Q(5.0 / 2.0, x2obs / 2.0);

        std::string random = "TRUE";
        if (pValue < 0.01) {
            random = "FALSE";
        }

        outputFile << v0 << ";&;" << v1 << ";&;" << v2 << ";&;" << v3 << ";&;" << v4 << ";&;" << v5 << ";&;" << x2obs
                   << ";&;" << pValue << ";&;" << random << ";\\\\" << std::endl;
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Operacja zakończona sukcesem." << std::endl;
}

int main() {
    FrequencyTestWithinABlock();
    FrequencyMonobitsTest();
    BinaryMatrixRankTest();
    RunsTest();
    LongestRunOfOnesInABlock();
    return 0;
}