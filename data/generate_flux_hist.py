import sys
import ROOT
import math

def gaus(x, width, x0):
    norm = 1. / math.sqrt(2*math.pi*width**2)
    return norm*math.exp(-(x-x0)**2 / (2*width**2))

def main(inputf, binw, output): 

    next(inputf)
    dat = [float(line.rstrip("\n")) for line in inputf]

    f = ROOT.TFile(output, "RECREATE")
    h = ROOT.TH1D("numu", "numu", len(dat), 0, len(dat)*binw)
    for i,d in enumerate(dat):
        # convert 1/m^2/10^6 POT -> 1/cm^2/10^20POT
        h.SetBinContent(i+1, d*1e-4*1e14)
    h.Write()
    f.Write()
    f.Close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("Usage: python generate_signal_response.py <INPUT.txt> <BINWIDTH>")
    else:
        with open(sys.argv[1]) as f:
            main(f, float(sys.argv[2]), "flux.root")
