"""
Script to extract rho00 parameter from simulated D* cost* distributions 
"""

import argparse
import numpy as np
import uproot
from ROOT import TF1, TFile, TH1F, gROOT, gStyle, TGaxis, TLine, TCanvas, TLegend, kAzure, kRed, kBlack, kGray

gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.035)
gStyle.SetTitleOffset(1.5, 'y')
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description="Arguments")
parser.add_argument("inputFile", metavar="text", default="AnalysisResults.root")
parser.add_argument("outputFile", metavar="text", default="DstarPolarization.root")
parser.add_argument("--batch", action="store_true", default=False)
args = parser.parse_args()

if args.batch:
    gROOT.SetBatch(True)

ptMins = [0, 1, 2, 3, 4, 5, 7, 10]
ptMaxs = [1, 2, 3, 4, 5, 7, 10, 20]
ptLims = ptMins.copy()
ptLims.append(ptMaxs[-1])
ptLims = np.array(ptLims, "f")

inFile = TFile.Open(args.inputFile)
treeDstar = uproot.open(args.inputFile)["tupleDstars"]
dfDstar = treeDstar.arrays(library="pd")
dfDstarNonPrompt = dfDstar.query("fromB == 1")
dfDstarPrompt = dfDstar.query("fromB == 0")

hRho00Prompt = TH1F("hRho00Prompt", ";#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}", len(ptMins), ptLims)
hRho00Prompt.SetLineWidth(2)
hRho00Prompt.SetLineColor(kRed+1)
hRho00Prompt.SetMarkerStyle(20)
hRho00Prompt.SetMarkerColor(kRed+1)
hRho00NonPrompt = TH1F("hRho00NonPrompt", ";#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}", len(ptMins), ptLims)
hRho00NonPrompt.SetLineWidth(2)
hRho00NonPrompt.SetLineColor(kAzure+4)
hRho00NonPrompt.SetMarkerStyle(20)
hRho00NonPrompt.SetMarkerColor(kAzure+4)

hCosThetaStarPrompt, hCosThetaStarNonPrompt, fCosThetaStarPrompt, fCosThetaStarNonPrompt = ([] for _ in range(4))
outFile = TFile(args.outputFile, "recreate")
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    hCosThetaStarPrompt.append(
        TH1F(
            f"hCosThetaStarPrompt_pt{ptMin}-{ptMax}",
            ";cos(#theta*);d#it{N}/dcos(#theta*)",
            10,
            0.,
            1.
        )
    )
    hCosThetaStarNonPrompt.append(
        TH1F(
            f"hCosThetaStarNonPrompt_pt{ptMin}-{ptMax}",
            ";cos(#theta*);d#it{N}/dcos(#theta*)",
            10,
            0.,
            1.
        )
    )

    fCosThetaStarPrompt.append(
        TF1(
            f"fCosThetaStarPrompt_pt{ptMin}-{ptMax}",
            "[0] * ( (1-[1]) + (3*[1] - 1) * x*x)",
            0.,
            1.
        )
    )
    fCosThetaStarNonPrompt.append(
        TF1(
            f"fCosThetaStarNonPrompt_pt{ptMin}-{ptMax}",
            "[0] * ( (1-[1]) + (3*[1] - 1) * x*x)",
            0.,
            1.
        )
    )
    hCosThetaStarPrompt[iPt].SetLineWidth(2)
    hCosThetaStarPrompt[iPt].SetLineColor(kBlack)
    hCosThetaStarPrompt[iPt].SetMarkerStyle(20)
    hCosThetaStarPrompt[iPt].SetMarkerColor(kBlack)

    hCosThetaStarNonPrompt[iPt].SetLineWidth(2)
    hCosThetaStarNonPrompt[iPt].SetLineColor(kBlack)
    hCosThetaStarNonPrompt[iPt].SetMarkerStyle(20)
    hCosThetaStarNonPrompt[iPt].SetMarkerColor(kBlack)

    fCosThetaStarPrompt[iPt].SetLineWidth(2)
    fCosThetaStarPrompt[iPt].SetLineColor(kAzure+4)

    fCosThetaStarNonPrompt[iPt].SetLineWidth(2)
    fCosThetaStarNonPrompt[iPt].SetLineColor(kAzure+4)

    dfDstarPromptPt = dfDstarPrompt.query(f"{ptMin} < pTDstar < {ptMax}")
    for cost in dfDstarPromptPt["cosThetaStar"].to_numpy():
        hCosThetaStarPrompt[iPt].Fill(cost)

    hCosThetaStarPrompt[iPt].Fit(fCosThetaStarPrompt[iPt], "Q")
    hRho00Prompt.SetBinContent(iPt+1, fCosThetaStarPrompt[iPt].GetParameter(1))
    hRho00Prompt.SetBinError(iPt+1, fCosThetaStarPrompt[iPt].GetParError(1))

    dfDstarNonPromptPt = dfDstarNonPrompt.query(f"{ptMin} < pTDstar < {ptMax}")
    for cost in dfDstarNonPromptPt["cosThetaStar"].to_numpy():
        hCosThetaStarNonPrompt[iPt].Fill(cost)

    hCosThetaStarNonPrompt[iPt].Fit(fCosThetaStarNonPrompt[iPt], "Q")
    hRho00NonPrompt.SetBinContent(iPt+1, fCosThetaStarNonPrompt[iPt].GetParameter(1))
    hRho00NonPrompt.SetBinError(iPt+1, fCosThetaStarNonPrompt[iPt].GetParError(1))

    outFile.cd()
    hCosThetaStarPrompt[iPt].Write()
    hCosThetaStarNonPrompt[iPt].Write()

lineRho00 = TLine(ptMins[0], 1./3, ptMaxs[-1], 1./3)
lineRho00.SetLineWidth(2)
lineRho00.SetLineColor(kGray+2)
lineRho00.SetLineStyle(9)

cRho00 = TCanvas("cRho00", "", 500, 500)
legRho00 = TLegend(0.2, 0.2, 0.4, 0.35)
legRho00.SetTextSize(0.045)
legRho00.SetBorderSize(0)
legRho00.SetFillStyle(0)
legRho00.SetHeader("PYTHIA8 + EvtGen")
legRho00.AddEntry(hRho00Prompt, "c #rightarrow D*^{+}", "p")
legRho00.AddEntry(hRho00NonPrompt, "b #rightarrow D*^{+}", "p")

hFrame = cRho00.DrawFrame(
    ptMins[0],
    0.2,
    ptMaxs[-1],
    0.5,
    ";#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}"
)
lineRho00.Draw("same")
hFrame.GetYaxis().SetDecimals()
hRho00Prompt.DrawCopy("same")
hRho00NonPrompt.DrawCopy("same")
legRho00.Draw()

outFile.cd()
hRho00Prompt.Write()
hRho00NonPrompt.Write()
cRho00.Write()
outFile.Close()

cRho00.SaveAs(args.outputFile.replace(".root", ".pdf"))

if not args.batch:
    input("Press enter to exit")
