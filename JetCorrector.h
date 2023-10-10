#ifndef JETCORRECTOR_H
#define JETCORRECTOR_H
// SingleJetCorrector
// v3.0
// Author: Yi Chen
// 
// This class applies JEC for any given level using TF1 as the workhorse
// Supposedly runs faster than v1.0
// v3.0: one can add list of text files to apply them one by one

// С++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// ROOT headers
#include "TObject.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

//________________
class SingleJetCorrector {
  public:
    SingleJetCorrector();
    SingleJetCorrector(std::string File);
    virtual ~SingleJetCorrector();
    enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
    void SetJetPT(double value)     { JetPT = value; }
    void SetJetEta(double value)    { JetEta = value; }
    void SetJetPhi(double value)    { JetPhi = value; }
    void SetJetArea(double value)   { JetArea = value; }
    void SetRho(double value)       { Rho = value; }
    void Initialize(std::string FileName);
    std::vector<std::string> BreakIntoParts(std::string Line);
    bool CheckDefinition(std::string Line);
    std::string StripBracket(std::string Line);
    SingleJetCorrector::Type ToType(std::string Line);
    double GetCorrection();
    double GetCorrectedPT();
    double GetValue(Type T);
  private:
    std::string Hack4(std::string Formula, char V, int N);
    bool Initialized;
    bool IsFunction;
    double JetPT, JetEta, JetPhi, JetArea, Rho;
    std::vector<std::string> Formulas;
    std::vector<std::vector<double>> Parameters;
    std::vector<std::vector<Type>> BinTypes;
    std::vector<std::vector<double>> BinRanges;
    std::vector<std::vector<Type>> Dependencies;
    std::vector<std::vector<double>> DependencyRanges;
    std::vector<TF1 *> Functions;

    ClassDef(SingleJetCorrector, 0)
};

//________________
class JetCorrector {
  public:
    JetCorrector()                               {}
    JetCorrector(std::string File)               { Initialize(File); }
    JetCorrector(std::vector<std::string> Files) { Initialize(Files); }
    void Initialize(std::string File)            { std::vector<std::string> X; X.push_back(File); Initialize(X); }
    void Initialize(std::vector<std::string> Files);
    void SetJetPT(double value)     { JetPT = value; }
    void SetJetEta(double value)    { JetEta = value; }
    void SetJetPhi(double value)    { JetPhi = value; }
    void SetJetArea(double value)   { JetArea = value; }
    void SetRho(double value)       { Rho = value; }
    double GetCorrection();
    double GetCorrectedPT();
  private:
    std::vector<SingleJetCorrector> JEC;
    double JetPT, JetEta, JetPhi, JetArea, Rho;

    ClassDef(JetCorrector, 0)
};

#endif // #define JETCORRECTOR_H