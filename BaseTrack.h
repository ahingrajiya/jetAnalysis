/**
 * @file BaseTrack.h
 * @author Abhishek Hingrajiya (abhi.h123@gmail.com)
 * @brief Base class for track descriptions
 * @version 0.1
 * @date 2024-06-28
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BaseTrack_h
#define BaseTrack_h

// ROOT Headers

#include "TObject.h"
#include "TVector3.h"

class BaseTrack : public TObject
{
public:
    ///@brief Constructor
    BaseTrack();
    ///@brief Destructor
    virtual ~BaseTrack()
    {
    }
    ///
    /// Setters
    ///
    ///@brief  Set track transverse momentum
    void setTrkPt(const Float_t &pt) { fTrkPt = {pt}; }
    ///@brief Set track pseudorapidity
    void setTrkEta(const Float_t &eta) { fTrkEta = {eta}; }
    ///@brief Set track phi
    void setTrkPhi(const Float_t &phi) { fTrkPhi = {phi}; }
    ///@brief Set Track Charge
    void setTrkChg(const Int_t &chg) { fTrkChg = {chg}; }
    ///
    /// Getter
    ///
    ///@brief Transverse momentum
    Float_t TrkPt() const { return fTrkPt; }
    ///@brief Pseudorapidity
    Float_t TrkEta() const { return fTrkEta; }
    ///@brief Azimuthal Angle
    Float_t TrkPhi() const { return fTrkPhi; }
    ///@brief Return three vectors of tracks
    TVector3 vec() const
    {
        TVector3 v;
        v.SetPtEtaPhi(fTrkPt, fTrkEta, fTrkPhi);
        return v;
    }

private:
    ///@brief  Track  Transverse Momentum
    Float_t fTrkPt;
    ///@brief Track Pseudorapidity
    Float_t fTrkEta;
    ///@brief Track azimuthal angle
    Float_t fTrkPhi;
    ///@brief Track Charge
    Int_t fTrkChg;

    ClassDef(BaseTrack, 1)
};
#endif