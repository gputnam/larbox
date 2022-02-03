////////////////////////////////////////////////////////////////////////
// Class:       Digitization
// Plugin Type: producer (Unknown Unknown)
// File:        Digitization_module.cc
//
// Generated at Mon Jan 31 14:30:07 2022 by Gray Putnam using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Data Objects
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Wire.h"

// LArSoft stuff
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "lardata/Utilities/AssociationUtil.h"

// Random numbers
#include "CLHEP/Random/RandFlat.h"

// ROOT
#include "TH1D.h"
#include "TFile.h"
#include "TComplex.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"

#include <memory>

namespace larbox {
  class Digitization;
}


class larbox::Digitization : public art::EDProducer {
public:
  explicit Digitization(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Digitization(Digitization const&) = delete;
  Digitization(Digitization&&) = delete;
  Digitization& operator=(Digitization const&) = delete;
  Digitization& operator=(Digitization&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Internal Functions
  void MakeNoise(std::vector<double> &v);
  void MakeSignal(const std::vector<double> &s, std::vector<double> &v);

  ~Digitization();

private:

  // Configuration

  // Input Data
  art::InputTag fSimChannelLabel;

  // Noise Configuration
  std::string fNoiseFileName;
  std::string fNoiseHistoName;
  std::vector<double> fNoiseSpectrum;
  double fNoiseFreqRandomness;
  double fNoiseScale;

  // Signal Configuration
  std::string fSignalFileName;
  std::string fSignalHistoName;
  std::vector<double> fSignalShape;
  double fSignalGain;

  // Data output Configuration
  bool fProduceDigits;
  bool fProduceWires;
  int fWireROIRange;

  // Random numbers
  CLHEP::HepRandomEngine &fEngine;

  // FFT algorithms
  TFFTRealComplex *fFFT;
  TFFTComplexReal *fInvFFT;

};


larbox::Digitization::Digitization(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fSimChannelLabel(p.get<art::InputTag>("SimChannelLabel")),
    fNoiseFileName(p.get<std::string>("NoiseFileName")),
    fNoiseHistoName(p.get<std::string>("NoiseHistoName")),
    fNoiseFreqRandomness(p.get<double>("NoiseFreqRandomness")),
    fNoiseScale(p.get<double>("NoiseScale")),
    fSignalFileName(p.get<std::string>("SignalFileName")),
    fSignalHistoName(p.get<std::string>("SignalHistoName")),
    fSignalGain(p.get<double>("SignalGain")),
    fProduceDigits(p.get<bool>("ProduceDigits")),
    fProduceWires(p.get<bool>("ProduceWires")),
    fWireROIRange(p.get<int>("WireROIRange")),
    fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "digitization", p, "Seed"))
{
  // Load the Noise Histogram
  {
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fNoiseFileName, fullFileName);

    TFile inputFile(fullFileName.c_str(), "READ");
    TH1D* noiseHist = (TH1D*)inputFile.Get(fNoiseHistoName.c_str());

    for (int ibin = 0; ibin < noiseHist->GetNbinsX(); ibin++) {
      fNoiseSpectrum.push_back(noiseHist->GetBinContent(ibin+1));
    }
  }

  // Load the Signal Response
  {
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fSignalFileName, fullFileName);

    TFile inputFile(fullFileName.c_str(), "READ");
    TH1D* signalHist = (TH1D*)inputFile.Get(fSignalHistoName.c_str());

    for (int ibin = 0; ibin < signalHist->GetNbinsX(); ibin++) {
      fSignalShape.push_back(signalHist->GetBinContent(ibin+1));
    }
  }

  // Setup FFT's
  auto const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  unsigned n_time_samples = detprop.NumberTimeSamples();
  fFFT = new TFFTRealComplex(n_time_samples, false);
  fFFT->Init("M", 0, NULL);
  fInvFFT = new TFFTComplexReal(n_time_samples, false);
  fInvFFT->Init("M", 0, NULL);

  // Output Info
  if (fProduceDigits) {
    produces<std::vector<raw::RawDigit>>();
  }

  if (fProduceWires) {
    produces<std::vector<recob::Wire>>();
  }

  if (fProduceDigits && fProduceWires) {
    produces<art::Assns<raw::RawDigit, recob::Wire>>();
  }

}

larbox::Digitization::~Digitization() {
  delete fFFT;
  delete fInvFFT;
}

void larbox::Digitization::produce(art::Event& e)
{
  // Output Data
  std::unique_ptr<std::vector<raw::RawDigit>> digits(new std::vector<raw::RawDigit>);
  std::unique_ptr<std::vector<recob::Wire>> wires(new std::vector<recob::Wire>);
  std::unique_ptr<art::Assns<raw::RawDigit, recob::Wire>> assn(new art::Assns<raw::RawDigit, recob::Wire>);

  art::PtrMaker<raw::RawDigit> DigitPtrMaker {e};
  art::PtrMaker<recob::Wire> WirePtrMaker {e};

  // Services
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>(); 
  auto const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  unsigned n_time_samples = detprop.NumberTimeSamples();

  // Input Data
  art::ValidHandle<std::vector<sim::SimChannel>> simChannels = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);

  std::map<unsigned, const sim::SimChannel *> simChannel_map;
  for (const sim::SimChannel &sc: *simChannels) {
    simChannel_map[sc.Channel()] = &sc;
  }

  // Get all the channels in the detector
  std::vector<raw::ChannelID_t> allChannels = geo->ChannelsInTPCs();

  for (unsigned c: allChannels) {
    std::vector<double> noise(n_time_samples, 0.);
    std::vector<double> signal(n_time_samples, 0.);
    std::vector<double> charge(n_time_samples, 0.);

    std::vector<int> charge_ticks;

    if (simChannel_map.count(c)) { // TODO: suppress channels with no charge?
      const sim::SimChannel &sc = *simChannel_map.at(c);
      for (unsigned tick = 0; tick < n_time_samples; tick++) {
        int tdc = clockData.TPCTick2TDC(tick);     

        charge[tick] += sc.Charge(tdc) / fSignalGain;

        if (sc.Charge(tdc) > 1e-4) {
          charge_ticks.push_back(tick);
        }
      }

      MakeSignal(charge, signal);
    }

    MakeNoise(noise);

    std::vector<short> adcs(n_time_samples, 0);

    for (unsigned tick = 0; tick < n_time_samples; tick++) {
      adcs[tick] = std::round(noise[tick] + signal[tick]);
    }

    raw::RawDigit d(c, n_time_samples, adcs);

    if (fProduceDigits) {
      digits->push_back(d);
    }

    // Make the wire if configured
    if (fProduceWires) {
      recob::Wire::RegionsOfInterest_t roiVec;
      
      for (int tick: charge_ticks) {
        std::vector<short> thisrange;
        std::copy(adcs.begin() + std::max((int)0, tick - fWireROIRange), adcs.begin() + std::min((int)adcs.size(), tick + fWireROIRange), 
            std::back_inserter(thisrange));

        roiVec.add_range(std::max((int)0, tick - fWireROIRange), std::move(thisrange));
      }

      recob::Wire w(std::move(roiVec), c, geo->View(c));

      wires->push_back(w);

      // Make the association if configured
      if (fProduceDigits) {
        art::Ptr<raw::RawDigit> ptr_digit = DigitPtrMaker(digits->size() - 1);
        art::Ptr<recob::Wire> ptr_wire = WirePtrMaker(wires->size() - 1);
        assn->addSingle(ptr_digit, ptr_wire);
      }
    }

  }

  if (fProduceDigits) {
    e.put(std::move(digits));
  }

  if (fProduceWires) {
    e.put(std::move(wires));
  }

  if (fProduceDigits && fProduceWires) {
    e.put(std::move(assn));
  }

}

void larbox::Digitization::MakeSignal(const std::vector<double> &s, std::vector<double> &v) {
  unsigned n_freq_points = s.size()/2 + 1;

  // Transform the Signal
  for (unsigned i = 0; i < s.size(); i++) {
    fFFT->SetPoint(i, s[i]);
  }
  fFFT->Transform();

  std::vector<TComplex> s_freq;
  for (unsigned i = 0; i < n_freq_points; i++) {
    double re, im;
    fFFT->GetPointComplex(i, re, im);
    s_freq.push_back(TComplex(re, im));
  }

  // And the Kernel
  for (unsigned i = 0; i < fSignalShape.size(); i++) {
    fFFT->SetPoint(i, fSignalShape[i]);
  }
  fFFT->Transform();

  std::vector<TComplex> k_freq;
  for (unsigned i = 0; i < n_freq_points; i++) {
    double re, im;
    fFFT->GetPointComplex(i, re, im);
    k_freq.push_back(TComplex(re, im));
  }

  // Multiply
  std::vector<TComplex> freq;
  for (unsigned i = 0; i < s_freq.size(); i++) {
    freq.push_back(s_freq[i]*k_freq[i]);
  }

  // Transform back to Time domain
  for (unsigned i = 0; i < freq.size(); i++) {
    fInvFFT->SetPointComplex(i, freq[i]);
  }
  fInvFFT->Transform();

  // Save to output
  for (unsigned i = 0; i < v.size(); i++) {
    v[i] = fInvFFT->GetPointReal(i, false) / v.size();
  }
}
 

void larbox::Digitization::MakeNoise(std::vector<double> &v) {
  // Load-up Random Numbers
  CLHEP::RandFlat flat(fEngine,-1,1);

  // Get the size of the noise in freq. space
  int n_freq_points = v.size()/2 + 1;

  std::vector<TComplex> noisef(n_freq_points, 0.);

  for (unsigned i = 0; i < noisef.size(); i++) {
    double r = flat.shoot();
    double thisf = fNoiseSpectrum.at(i) * (1 - fNoiseFreqRandomness + 2*fNoiseFreqRandomness * r) * fNoiseScale; 
    double phase = flat.shoot() * 2 * M_PI;

    TComplex fc(thisf * cos(phase), thisf * sin(phase));
    noisef[i] = fc;
  }

  // Fourier Transform to the output
  for (unsigned i = 0; i < noisef.size(); i++) {
    fInvFFT->SetPointComplex(i, noisef[i]);
  }
  fInvFFT->Transform();

  // Get output
  for (unsigned i = 0; i < v.size(); i++) {
    v[i] = fInvFFT->GetPointReal(i, false);
  }  

}

DEFINE_ART_MODULE(larbox::Digitization)
