/**
 * @file GenTrack.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Generated track description
 * @version 0.1
 * @date 2023-10-19
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef GenTrack_h
#define GenTrack_h

// ROOT headers
#include "TObject.h"
#include "BaseTrack.h"

//_________________
class GenTrack : public BaseTrack
{
public:
  /// @brief constructor
  GenTrack();
  /// @brief destructor
  virtual ~GenTrack() {}

  /// Setter

  ///@brief Set Particle Identification
  void setPDGID(const Int_t &pdg) { fPDGID = {pdg}; }
  ///@brief Set Subevent
  void setSube(const Char_t &sube) { fSube = {sube}; }
  ///@brief Print parameters of tracks
  void print();

  // Getter

  ///@brief Get PDGID
  Int_t TrkPDGID() const { return fPDGID; }
  ///@brief Get Sube
  Char_t TrkSube() const { return fSube; }

private:
  ///@brief Particle ID
  Int_t fPDGID;
  ///@brief Subevent Identification
  Char_t fSube;

  ClassDef(GenTrack, 1)
};

#endif // #define GenTrack_h