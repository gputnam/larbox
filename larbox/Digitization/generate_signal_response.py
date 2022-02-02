import sys
import ROOT
import math

def gaus(x, width, x0):
    norm = 1. / math.sqrt(2*math.pi*width**2)
    return norm*math.exp(-(x-x0)**2 / (2*width**2))

def main(nticks, width, output):
    f = ROOT.TFile(output, "RECREATE")
    h = ROOT.TH1D("Response", "Response", nticks, 0, nticks)
    for i in range(nticks):
        h.SetBinContent(i+1, gaus(min(i, nticks-i), width, 0))
    h.Write()
    f.Write()
    f.Close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("Usage: python generate_signal_response.py <NTICK> <WIDTH>")
    else:
        main(int(sys.argv[1]), float(sys.argv[2]), "response.root")
