"""
Script to extract rho00 parameter from simulated D* cost* distributions 
"""

import argparse
import numpy as np
import uproot
from ROOT import TF1, TFile, TH1F, gROOT, gStyle, TGaxis, TLine, TCanvas, TLegend
from ROOT import TF1, TFile, TH1F, gROOT, gStyle, TGaxis, TLine, TCanvas, TLegend, kAzure, kBlue, kRed, kBlack, kGray, kRainBow

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
gStyle.SetPalette(kRainBow)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description="Arguments")
parser.add_argument("inputFile", metavar="text", default="AnalysisResults.root")
parser.add_argument("outputFile", metavar="text", default="DstarPolarization.root")
parser.add_argument("--batch", action="store_true", default=False)
parser.add_argument("--ispythia", action="store_true", default=False)
args = parser.parse_args()

if args.batch:
    gROOT.SetBatch(True)

ptMins = [0, 1, 2, 3, 4, 5, 7, 10]
ptMaxs = [1, 2, 3, 4, 5, 7, 10, 20]
ptLims = ptMins.copy()
ptLims.append(ptMaxs[-1])
ptLims = np.array(ptLims, "f")

inFile = TFile.Open(args.inputFile)
hPxCorr = inFile.Get("hPxCorr")
hPyCorr = inFile.Get("hPyCorr")
hPzCorr = inFile.Get("hPzCorr")
hEtaCorr = inFile.Get("hEtaCorr")
hPhiCorr = inFile.Get("hPhiCorr")
hPxCorr.SetDirectory(0)
hPyCorr.SetDirectory(0)
hPzCorr.SetDirectory(0)
hEtaCorr.SetDirectory(0)
hPhiCorr.SetDirectory(0)
inFile.Close()

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
hRho00BNonPrompt = TH1F("hRho00BNonPrompt", ";#it{p}_{T} (GeV/#it{c});#it{#rho}_{00} (B frame)", len(ptMins), ptLims)
hRho00BNonPrompt.SetLineWidth(2)
hRho00BNonPrompt.SetLineColor(kBlue+3)
hRho00BNonPrompt.SetMarkerStyle(20)
hRho00BNonPrompt.SetMarkerColor(kBlue+3)

hCosThetaStarPrompt, hCosThetaStarNonPrompt, hCosThetaStarBNonPrompt, \
    fCosThetaStarPrompt, fCosThetaStarNonPrompt, fCosThetaStarBNonPrompt = ([] for _ in range(6))
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
    hCosThetaStarBNonPrompt.append(
        TH1F(
            f"hCosThetaStarBNonPrompt_pt{ptMin}-{ptMax}",
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
    fCosThetaStarBNonPrompt.append(
        TF1(
            f"fCosThetaStarBNonPrompt_pt{ptMin}-{ptMax}",
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

    hCosThetaStarBNonPrompt[iPt].SetLineWidth(2)
    hCosThetaStarBNonPrompt[iPt].SetLineColor(kBlack)
    hCosThetaStarBNonPrompt[iPt].SetMarkerStyle(20)
    hCosThetaStarBNonPrompt[iPt].SetMarkerColor(kBlack)

    fCosThetaStarPrompt[iPt].SetLineWidth(2)
    fCosThetaStarPrompt[iPt].SetLineColor(kAzure+4)

    fCosThetaStarNonPrompt[iPt].SetLineWidth(2)
    fCosThetaStarNonPrompt[iPt].SetLineColor(kAzure+4)

    fCosThetaStarBNonPrompt[iPt].SetLineWidth(2)
    fCosThetaStarBNonPrompt[iPt].SetLineColor(kAzure+4)

    dfDstarPromptPt = dfDstarPrompt.query(f"{ptMin} < pTDstar < {ptMax}")
    for cost in dfDstarPromptPt["cosThetaStar"].to_numpy():
        hCosThetaStarPrompt[iPt].Fill(cost)

    hCosThetaStarPrompt[iPt].Fit(fCosThetaStarPrompt[iPt], "Q")
    hRho00Prompt.SetBinContent(iPt+1, fCosThetaStarPrompt[iPt].GetParameter(1))
    hRho00Prompt.SetBinError(iPt+1, fCosThetaStarPrompt[iPt].GetParError(1))

    dfDstarNonPromptPt = dfDstarNonPrompt.query(f"{ptMin} < pTDstar < {ptMax}")
    for cost in dfDstarNonPromptPt['cosThetaStar'].to_numpy():
        hCosThetaStarNonPrompt[iPt].Fill(cost)

    hCosThetaStarNonPrompt[iPt].Fit(fCosThetaStarNonPrompt[iPt], "Q")
    hRho00NonPrompt.SetBinContent(iPt+1, fCosThetaStarNonPrompt[iPt].GetParameter(1))
    hRho00NonPrompt.SetBinError(iPt+1, fCosThetaStarNonPrompt[iPt].GetParError(1))

    for costB in dfDstarNonPromptPt['cosThetaStarB'].to_numpy():
        hCosThetaStarBNonPrompt[iPt].Fill(costB)

    hCosThetaStarBNonPrompt[iPt].Fit(fCosThetaStarBNonPrompt[iPt], "Q")
    hRho00BNonPrompt.SetBinContent(iPt+1, fCosThetaStarBNonPrompt[iPt].GetParameter(1))
    hRho00BNonPrompt.SetBinError(iPt+1, fCosThetaStarBNonPrompt[iPt].GetParError(1))

    outFile.cd()
    hCosThetaStarPrompt[iPt].Write()
    hCosThetaStarNonPrompt[iPt].Write()
    hCosThetaStarBNonPrompt[iPt].Write()

lineRho00 = TLine(ptMins[0], 1./3, ptMaxs[-1], 1./3)
lineRho00.SetLineWidth(2)
lineRho00.SetLineColor(kGray+2)
lineRho00.SetLineStyle(9)

cRho00 = TCanvas("cRho00", "", 500, 500)
legRho00 = TLegend(0.2, 0.2, 0.4, 0.4)
legRho00.SetTextSize(0.045)
legRho00.SetBorderSize(0)
legRho00.SetFillStyle(0)
if args.ispythia:
    legRho00.SetHeader("PYTHIA8")
else:
    legRho00.SetHeader("PYTHIA8 + EvtGen")
legRho00.AddEntry(hRho00Prompt, "c #rightarrow D*^{+}", "p")
legRho00.AddEntry(hRho00NonPrompt, "b #rightarrow D*^{+}", "p")
legRho00.AddEntry(hRho00BNonPrompt, "b #rightarrow D*^{+} (B frame)", "p")

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
hRho00BNonPrompt.DrawCopy("same")
legRho00.Draw()

outFile.cd()
hRho00Prompt.Write()
hRho00NonPrompt.Write()
hRho00BNonPrompt.Write()
cRho00.Write()
outFile.Close()

cRho00.SaveAs(args.outputFile.replace(".root", ".pdf"))

cPCorr = TCanvas("cPCorr", "", 1200, 400)
cPCorr.Divide(3, 1)
cPCorr.cd(1)
hPxCorr.Draw("colz")
cPCorr.cd(2)
hPyCorr.Draw("colz")
cPCorr.cd(3)
hPzCorr.Draw("colz")

cAngleCorr = TCanvas("cAngleCorr", "", 800, 400)
cAngleCorr.Divide(2, 1)
cAngleCorr.cd(1)
hPhiCorr.Draw("colz")
cAngleCorr.cd(2)
hEtaCorr.Draw("colz")

cPCorr.SaveAs(args.outputFile.replace(".root", "_momCorr.pdf"))
cAngleCorr.SaveAs(args.outputFile.replace(".root", "_angleCorr.pdf"))

if not args.batch:
    input("Press enter to exit")
