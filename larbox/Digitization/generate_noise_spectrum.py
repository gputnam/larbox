import sys
import ROOT
import math

def main(nticks, output):
    f = ROOT.TFile(output, "RECREATE")
    h = ROOT.TH1D("NoiseFreq", "NoiseFreq", nticks, 0, nticks)
    for i in range(nticks):
        h.SetBinContent(i+1, 1)
    h.Write()
    f.Write()
    f.Close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print ("Usage: python generate_noise_spectrum.py <NTICK>")
    else:
        main(int(sys.argv[1])//2+1, "noise.root")
