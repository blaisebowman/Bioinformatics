#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <set>
#include <map>
#include <vector>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment 1

string multipleAlignments(int i, int j, vector<vector<string>> tbl, vector<vector<string>> p) {
//Function used to traceback and determine IF there is more than one optimal alignment
//Part F (all possible alignments is covered in the traceback function
// '\u2196' is the NORTHWEST arrow, '\u2191' is the UPWARDS arrow, '\u2190' is the LEFTWARDS arrow
    while ((i > 0 || j > 0) && p.size() <= 1) {
        if (tbl[i][j] == "\u2196") {
            i--, j--;
            //has not hit branch
        } else if (tbl[i][j] == "\u2191") {
            i--;
            //has not hit branch
        } else if (tbl[i][j] == "\u2190") {
            //has not hit branch
            j--;
        } else if (tbl[i][j] == "\u2196 \u2191") { //cell contains NW and UP arrows
            return "YES";
        } else if (tbl[i][j] == "\u2196 \u2190") { //cell contains NW and LEFT arrows
            return "YES";
        } else if (tbl[i][j] == "\u2191 \u2190") { //cell contains UP and LEFT arrows
            return "YES";
        } else if (tbl[i][j] == "\u2196 \u2191 \u2190") {//cell contains, NW, UP, and LEFT arrows
            return "YES";
        }
    }
    return "NO"; // has reached the beginning of the 2D vector without finding a cell with multiple path options
}


vector<vector<string>>
traceback(int i, int j, vector<vector<string>> tbl, string m, string n, string dir, vector<vector<string>> p, int v) {
//Function used to traceback and determine the optimal alignment(s), if there is more than one optimal alignment
// '\u2196' is the NORTHWEST arrow, '\u2191' is the UPWARDS arrow, '\u2190' is the LEFTWARDS arrow
// vector <vector <string>> p is the 2D vector used to keep track of the backtrace
// string dir is used to keep track of direction in the event that a cell contains multiple arrows
    while (i > 0 || j > 0) {
        if (tbl[i][j] == "\u2196" || dir == "\u2196") {
            p[v][0].insert(p[v][0].begin(), m[i - 1]); //inserts the i-1th char of Sequence 1 (n)
            p[v][1].insert(p[v][1].begin(), n[j - 1]); //inserts the i-1th char of Sequence 2 (m)
            dir = "";
            i--, j--; //move to the i-1th row and the j-1th column
        } else if (tbl[i][j] == "\u2191" || dir == "\u2191") {
            p[v][0].insert(p[v][0].begin(), m[i - 1]); //inserts the i-1th char of Sequence 1
            p[v][1].insert(p[v][1].begin(), '-'); //inserts a - into Sequence 2
            i--; //decrements the value of i, move to the i-1th row
            dir = "";
        } else if (tbl[i][j] == "\u2190" || dir == "\u2190") {
            p[v][0].insert(p[v][0].begin(), '-'); //inserts a - into Sequence 1
            p[v][1].insert(p[v][1].begin(), n[j - 1]); //the j-1th char of Sequence 2
            j--; //decrements the value of j, move to the previous j-1th column
            dir = "";
        } else if (tbl[i][j] == "\u2196 \u2191") { //cell contains NW and UP arrows
            vector<vector<string>> tmp = p; //copies the elements of p into tmp; needed to retain the contents of p
            p = traceback(i, j, tbl, m, n, "\u2196", tmp, v); //recursive call, evaluates path with the NW direction
            tmp.insert(tmp.end(), p.begin(), p.end()); //appends the result of the recursive call to 2D vector
            p = traceback(i, j, tbl, m, n, "\u2191", tmp, v); //recursive call, evaluates path with the UP direction
            return p;
        } else if (tbl[i][j] == "\u2196 \u2190") { //cell contains NW and LEFT arrows
            vector<vector<string>> tmp = p;
            p = traceback(i, j, tbl, m, n, "\u2196", tmp, v); //recursive call, evaluates path in the NW direction
            tmp.insert(tmp.end(), p.begin(), p.end());
            p = traceback(i, j, tbl, m, n, "\u2190", tmp, v); //recursive call, evaluates path in the LEFT direction
            return p;
        } else if (tbl[i][j] == "\u2191 \u2190") { //cell contains NW and LEFT arrows
            vector<vector<string>> tmp = p;
            p = traceback(i, j, tbl, m, n, "\u2191", tmp, v); //recursive call, evaluates path in the NW direction
            tmp.insert(tmp.end(), p.begin(), p.end());
            p = traceback(i, j, tbl, m, n, "\u2190", tmp, v); //recursive call, evaluates path in the LEFT direction
            return p;
        } else if (tbl[i][j] == "\u2196 \u2191 \u2190") {
            vector<vector<string>> tmp = p;
            p = traceback(i, j, tbl, m, n, "\u2196", tmp, v);
            tmp.insert(tmp.end(), p.begin(), p.end());
            p = traceback(i, j, tbl, m, n, "\u2191", tmp, v);
            tmp.insert(tmp.end(), p.begin(), p.end());
            p = traceback(i, j, tbl, m, n, "\u2190", tmp, v);
            return p;
        }
    }
    if (i == 0 && j == 0) {
        //if the we arrive at the first cell in the vector (upper-left corner), return p, the pair of sequence strings
        return p;
    }
}

int optimalScore(vector<vector<int>> const &matrix) {
    unsigned long row = matrix.size() - 1;
    unsigned long col = matrix[row].size() - 1;
    return matrix[row][col]; //returns the last entry in the 2D vector, which is the optimal score
}

pair<vector<vector<int>>, vector<vector<string>>> vectPair(string const &m, const string &n) {
    vector<vector<int>> matrix(m.length() + 1, vector<int>(n.length() + 1, 0));
    vector<vector<string>> aws(m.length() + 1, vector<string>(n.length() + 1, "*"));
    for (int i = 0; i < matrix.size(); i++) { //row
        matrix[i][0] = i * -2;
        for (int j = 0; j < matrix[i].size(); j++) {
            matrix[0][j] = j * -2;
        }
    }
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (i == 0 && j == 0) {
                aws[i][j] = "*";
            } else if (i == 0 && j >= 0) {
                aws[i][j] = "\u2190";
            } else if (j == 0 && i >= 0) {
                aws[i][j] = "\u2191";
            } else {
                int val = (m[i - 1] == n[j - 1]) ? 2 : -1;
                matrix[i][j] = max(((matrix[i - 1][j - 1] + val)),
                                   max(((matrix[i - 1][j] - 2)), ((matrix[i][j - 1] - 2))));
                if (matrix[i][j] == matrix[i - 1][j - 1] + val) {
                    if ((i < m.length() + 1 && j < n.length() + 1) &&
                        (matrix[i - 1][j - 1] + val > matrix[i][j - 1] - 2) &&
                        (matrix[i - 1][j - 1] + val > matrix[i - 1][j] - 2)) {
                        aws[i][j] = "\u2196"; //DIAGONAL
                    } else if ((i < m.length() + 1 && j < n.length() + 1) &&
                               (matrix[i - 1][j - 1] + val == matrix[i][j - 1] - 2) &&
                               (matrix[i - 1][j - 1] + val > matrix[i - 1][j] - 2)) {
                        aws[i][j] = "\u2196 \u2190"; //Multiple paths available
                    } else if ((i < m.length() + 1 && j < n.length() + 1) &&
                               (matrix[i - 1][j - 1] + val != matrix[i][j - 1] - 2) &&
                               (matrix[i - 1][j - 1] + val == matrix[i - 1][j] - 2)) {
                        aws[i][j] = "\u2196 \u2191"; //Multiple paths available
                    } else {
                        aws[i][j] = "\u2196 \u2191 \u2190"; //Multiple paths available
                    }

                } else if (matrix[i][j] == matrix[i - 1][j] - 2) {
                    if ((i < m.length() + 1 && j < n.length() + 1) &&
                        (matrix[i - 1][j] - 2 > matrix[i - 1][j - 1] + val) &&
                        (matrix[i - 1][j] - 2 > matrix[i][j - 1] - 2)) {
                        aws[i][j] = "\u2191";
                    } else if ((i < m.length() + 1 && j < n.length() + 1) &&
                               (matrix[i - 1][j] - 2 == matrix[i - 1][j - 1] + val) &&
                               (matrix[i - 1][j] - 2 > matrix[i][j - 1] - 2)) {
                        aws[i][j] = "\u2196 \u2191"; //Multiple paths available
                    } else {
                        aws[i][j] = "\u2191 \u2190"; //Multiple paths available
                    }
                } else if (matrix[i][j] == matrix[i][j - 1] - 2) {
                    if ((i < m.length() + 1 && j < n.length() + 1) &&
                        (matrix[i][j - 1] - 2 > matrix[i - 1][j - 1] + val) &&
                        (matrix[i][j - 1] - 2 > matrix[i - 1][j] - 2)) {
                        aws[i][j] = "\u2190";
                    } else if ((i <= m.length() + 1 && j <= n.length() + 1) &&
                               (matrix[i][j - 1] - 2 == matrix[i - 1][j - 1] + val) &&
                               (matrix[i - 1][j] - 2 != matrix[i][j - 1] - 2)) {
                        aws[i][j] = "\u2196 \u2190"; //Multiple paths available
                    } else {
                        aws[i][j] = "\u2191 \u2190"; //Multiple paths available
                    }
                }
            }
        }
    }
    return make_pair(matrix, aws);
}

pair<string, string> defaultAlignment(vector<vector<int>> const &dp, string m, string n) {
    int i = m.length();
    int j = n.length();
    string a, b;
    while (i > 0 && j > 0) {
        int val = (m[i - 1] == n[j - 1]) ? 2 : -1;
        if ((i > 0) && (j > 0) && (dp[i][j] == dp[i - 1][j - 1] + val)) {
            a.insert(a.begin(), m[i - 1]);
            b.insert(b.begin(), n[j - 1]);
            i--;
            j--;
        } else if (i > 0 && (dp[i][j] == dp[i - 1][j] - 2)) {
            a.insert(a.begin(), m[i - 1]);
            b.insert(b.begin(), '-');
            i--;
        } else {
            a.insert(a.begin(), '-');
            b.insert(b.begin(), n[j - 1]);
            j--;
        }
    }
    return make_pair(a, b);
}

int main(int argc, char *argv[]) {
    // For testng within CLion 20.2.3, RUN -> EDIT CONFIGURATIONS -> PROGRAM ARGUMENTS
    // -> REDIRECT INPUT - path to .fasta file (assignment-1.fasta)
    //ifstream file ("assignment-1.fasta"); //for testing within CLion,

    ifstream file(argv[1]); // for testing via Linux / CMD
    ofstream outS;
    string line, line2;
    map<string, string> lines;
    vector<string> sequence;
    vector<string> data;
    if (file.is_open()) {
        string val;
        while (getline(file, line)) {
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.find('>') != string::npos) { //if the line begins with '>'
                sequence.push_back(line);
            } else {
                data.push_back(line);
            }
        }
        for (int i = 0; i < sequence.size(); i++) {
            lines.insert(pair<string, string>(sequence[i], data[i]));
        }
        file.close();
    }
    vector<string> sequences;
    for (const auto &s  : lines) {
        sequences.push_back(s.second);
    }
    string m = data[0]; //Sequence 1's data
    string n = data[1]; //Sequence 2's data

//Start Part A -> Write the optimal score as an integer
    ofstream ofile("assignment1.o1");
    pair<vector<vector<int>>, vector<vector<string>>> matrices = vectPair(m, n);
    if (ofile.is_open()) {
        ofile << optimalScore(matrices.first); //writes the optimal score to assignment1.o1
    }
    ofile.close();
    ofile.clear();
//END Part A

//Start Part B -> Display the dynamic programming matrix for the pair of sequences
    ofile.open("assignment1.o2");
    if (ofile.is_open()) {
        for (const auto &row : matrices.first) { //writes matrix to assignment1.o2
            for (const auto &s : row) {
                ofile << s << " ";
            }
            ofile << endl;
        }
    }
    ofile.close();
    ofile.clear();
//END Part B

//Start Part C -> Determine the best alignments of the pair of DNA sequences
    pair<string, string> alignment = defaultAlignment(matrices.first, m, n);
    ofile.open("assignment1.o3");
    if (ofile.is_open()) { //writes best alignment of the pair to to assignment1.o3
        ofile << alignment.first << endl;
        ofile << alignment.second << endl;
    }
    ofile.close();
    ofile.clear();
//END Part C


//Start Part D -> Is there more than one optimal alignment?
    vector<vector<string>> multiple;
    multiple.resize(1, vector<string>(2, ""));
    string alignments = multipleAlignments(m.length(), n.length(), matrices.second, multiple);
    ofile.open("assignment1.o4");
    if (ofile.is_open()) { //writes YES or NO to assignment1.o4
        if (alignments == "YES") {
            ofile << "YES" << endl;
        } else {
            ofile << "NO" << endl;
        }
    }
    ofile.close();
    ofile.clear();
//END Part D

//Start Part F -> Report ALL possible alignments.

/* To SAFELY perform algorithm analysis (Part E) on determining ALL optimal alignments (Part E) for longer pairs of DNA sequences,
remove the // preceding the /* above vector <vector <string>> p and remove the // and above the line return 0, remove the preceding the */

/*
Note 1: Algorithm Analysis did not mention finding all possible alignments.
Note 2: Algorithm Analysis behaves as expected when determining the optimal alignment, and IF there is more than one optimal alignment.
This will report all optimal alignments, for pairs of DNA sequences up to 50 characters in length.
*/

//NOTE: for long strings, linux will kill the program
//Comment out this code if testing part E's sequences

    vector<vector<string>>
            p;
    p.resize(1, vector<string>(2, ""));
    vector<vector<string>> pairs = traceback(m.length(), n.length(), matrices.second, m, n,
                                             matrices.second[m.length() - 1][n.length() - 1], p, 0);
    ofile.open("assignment1.o5");
    if (ofile.is_open()) {
        for (const auto &row : pairs) { //writes ALL possible optimal alignments to assignment1.o5
            for (const auto &s : row) {
                ofile << s << "\n";
            }
            ofile << endl;
        }
    }
    ofile.close();
    ofile.clear();
    
    return 0;
    //END Part F
}