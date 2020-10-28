#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <set>
#include <random>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment 1 Sequence Generator

pair <string, string> algorithmAnalysis (const int sequenceLength){
    //given a sequence length, generate a pair of DNA sequences of length sequenceLength
    //length 50, 100, 250, 500, 1000, 1,500, 2000, etc.
    const string dnaChars ="ATCG";
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    uniform_int_distribution<>pos(0, dnaChars.size() - 1);
    string sequence1, sequence2;
    for (int i = 0; i < sequenceLength; i++){
        sequence1 += dnaChars[pos(randomGenerator)];
        sequence2 += dnaChars[pos(randomGenerator)];
    }
    if(sequence1 == sequence2){
        algorithmAnalysis(sequenceLength);
        //pairs of DNA sequences (NOT IDENTICAL)
    }
    return make_pair(sequence1, sequence2);
}

int main(int argc, char *argv []) {
    // Companion to bowman.cpp, used for Part E to generate random string sequence pairs

    ofstream ofile ("50characters.o1");
    pair <string, string> sequences50 = algorithmAnalysis(50);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences50.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences50.second << endl;
    }
    ofile.close();
    ofile.clear();

    ofile.open("100characters.o2");
    pair <string, string> sequences100 = algorithmAnalysis(100);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences100.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences100.second << endl;
    }
    ofile.close();
    ofile.clear();

    ofile.open("250characters.o3");
    pair <string, string> sequences250 = algorithmAnalysis(250);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences250.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences250.second << endl;
    }
    ofile.close();
    ofile.clear();

    ofile.open("500characters.o4");
    pair <string, string> sequences500 = algorithmAnalysis(500);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences500.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences500.second << endl;
    }
    ofile.close();
    ofile.clear();

    ofile.open("1000characters.o5");
    pair <string, string> sequences1000 = algorithmAnalysis(1000);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences1000.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences1000.second << endl;
    }
    ofile.close();
    ofile.clear();

    ofile.open("1500characters.o6");
    pair <string, string> sequences1500 = algorithmAnalysis(1500);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences1500.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences1500.second << endl;
    }
    ofile.close();
    ofile.clear();
    ofile.open("2000characters.o7");
    pair <string, string> sequences2000 = algorithmAnalysis(2000);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences2000.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences2000.second << endl;
    }
    ofile.close();
    ofile.clear();
    ofile.open("5000characters.o8");
    pair <string, string> sequences5000 = algorithmAnalysis(5000);
    if(ofile.is_open()) {
        ofile << "> Sequence 1" << endl;
        ofile << sequences5000.first << endl;
        ofile << "> Sequence 2" << endl;
        ofile << sequences5000.second << endl;
    }
    ofile.close();
    ofile.clear();
    return 0;
}

