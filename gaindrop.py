#!/usr/bin/env python
import math
import array
import ROOT
from ROOT import TFile, TH1D, TH2D, TCanvas, TPolyLine
import argparse
parser = argparse.ArgumentParser(description='gaindrop.py...')
parser.add_argument('-f', metavar='file',                                 required=True,  help='input file name')
parser.add_argument('-g', metavar='gap',                                  required=True,  help='gap number [1-4]')
parser.add_argument('-x', metavar='rectangle size in x [mm]',             required=False, help='rectangle size in x [mm]')
parser.add_argument('-y', metavar='rectangle size in y [mm]',             required=False, help='rectangle size in y [mm]')
parser.add_argument('-c', metavar='rectangle center x,y [mm]',            required=False, help='rectangle center x,y [mm]')
parser.add_argument('-z', metavar='current threshold out of the anomaly', required=False, help='current threshold out of the anomaly')
parser.add_argument('-l', metavar='log z?',                               required=False, help='log z?')
argus = parser.parse_args()
fname = argus.f
ngap  = argus.g

fnameout = fname
fnameout = fnameout.replace("root/","pdf/")
fnameout = fnameout.replace(".root","_")

Lx    = -1
Ly    = -1
x0    = -1
y0    = -1
zmin  = -1
logz  = False

### first iteration doesn't draw the rectangle
examine = (argus.x==None)

if(not examine):
   Lx    = float(argus.x)
   Ly    = float(argus.y)
   center = (argus.c).split(",")
   x0    = float(center[0])
   y0    = float(center[1])
   zmin  = float(argus.z)
   logz  = True if(argus.l=="1") else False

### rectangle corners
xmin = x0-Lx/2
xmax = x0+Lx/2
ymin = y0-Ly/2
ymax = y0+Ly/2

### some root stuff
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

### get the histos
f = TFile(fname,"READ")
h1X = f.Get("IvsX"+ngap)
h1Y = f.Get("IvsY"+ngap)
h2  = f.Get("IvsXY"+ngap)
h2.GetYaxis().UnZoom()
h2Anomal = h2.Clone("h2Anomal")
h2Anomal.Reset()
h2Full = h2.Clone("h2Full")
h2Full.Reset()

### draw the 
cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
if(logz): ROOT.gPad.SetLogz()
h2.Draw("colz")
h2.GetYaxis().UnZoom()
ROOT.gPad.UnZoomed()
cnv.SaveAs(fnameout+"gaindrop_2d.pdf")
if(examine):
   cnv.SaveAs(fnameout+"gaindrop.pdf")
   quit()
else:
   cnv.SaveAs(fnameout+"gaindrop.pdf(")

cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
h1X.Draw("hist")
cnv.SaveAs(fnameout+"gaindrop_1d.pdf")
cnv.SaveAs(fnameout+"gaindrop.pdf")

h1Xslice = h1X.Clone("sliceX")
h1Xslice.Reset()
by = h2.GetYaxis().FindBin(y0)
for bx in range(h1Xslice.GetNbinsX()+1):
   h1Xslice.SetBinContent(bx, h2.GetBinContent(bx,by))

cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
h1Xslice.Draw("hist")
cnv.SaveAs(fnameout+"gaindrop_slice.pdf")
cnv.SaveAs(fnameout+"gaindrop.pdf")



rectangle = TPolyLine( 5, array.array('d', [xmin,xmin,xmax,xmax,xmin]), array.array('d',[ymin,ymax,ymax,ymin,ymin]) )
rectangle.SetLineColor(ROOT.kWhite)

avg_rect = 0
pts_rect = 0
avg_rest = 0
pts_rest = 0
avg_blob = 0
pts_blob = 0
pts_full = 0
### first calculate the average in+out the initial rectangle
for bx in range(h2.GetXaxis().GetNbins()+1):
   x = h2.GetXaxis().GetBinCenter(bx)
   for by in range(h2.GetYaxis().GetNbins()+1):
      y = h2.GetYaxis().GetBinCenter(by)
      z = h2.GetBinContent(bx,by)
      if(x<xmin or x>xmax or y<ymin or y>ymax):
         if(z>zmin):
            avg_rest += z
            pts_rest += 1
      else: 
         avg_rect += z
         pts_rect += 1
avg_rect = avg_rect/pts_rect
avg_rest = avg_rest/pts_rest
drop_rect = (avg_rect-avg_rest)/avg_rest*100
### then get the contour area of the anomaly within the initial rectangle
for bx in range(h2.GetXaxis().GetNbins()+1):
   x = h2.GetXaxis().GetBinCenter(bx)
   for by in range(h2.GetYaxis().GetNbins()+1):
      y = h2.GetYaxis().GetBinCenter(by)
      z = h2.GetBinContent(bx,by)
      # if(z>0.05*avg_rest):
      if(z>zmin):
         pts_full += 1
         h2Full.Fill(x,y,50)
      if(x<xmin or x>xmax or y<ymin or y>ymax): continue
      if(z>=0.9*avg_rest): continue
      if(z<zmin): continue
      avg_blob += z
      pts_blob += 1
      h2Anomal.Fill(x,y,100)
avg_blob = avg_blob/pts_blob
drop_blob = (avg_blob-avg_rest)/avg_rest*100
area_blob = pts_blob*h2.GetXaxis().GetBinWidth(1)*h2.GetYaxis().GetBinWidth(1)
area_full = pts_full*h2.GetXaxis().GetBinWidth(1)*h2.GetYaxis().GetBinWidth(1)
# print("avg in problematic blob=%g, rect=%g, rest=%g --> drop(blob)=%g, drop(rect)=%g for blob area fraction=%.1f%%" % (avg_blob, avg_rect, avg_rest, drop_blob, drop_rect, area_blob/area_full*100))
#print("full area=%g [cm2], blob area=%g [cm2]" % (area_full,area_blob))

deepest = 1e10
for bx in range(h1Xslice.GetNbinsX()+1):
   x = h1Xslice.GetBinCenter(bx)
   if(x<xmin or x>xmax): continue
   z = h1Xslice.GetBinContent(bx)
   if(z>=avg_rest):      continue
   deepest = z if(z<deepest) else deepest
drop_singular = (deepest-avg_rest)/avg_rest*100
# print("avg in singular point=%g, rest=%g --> drop(singular)=%g%%" % (deepest, avg_rest, drop_singular))


cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
h2.Draw("colz")
h2.GetYaxis().UnZoom()
ROOT.gPad.UnZoomed()
if(logz): ROOT.gPad.SetLogz()
rectangle.Draw("same")
cnv.SaveAs(fnameout+"gaindrop_2d_anomaly.pdf")
cnv.SaveAs(fnameout+"gaindrop.pdf")

cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
h2.Draw("colz")
h2.GetYaxis().UnZoom()
ROOT.gPad.UnZoomed()
if(logz): ROOT.gPad.SetLogz()
h2Anomal.Draw("scat same")
h2Anomal.GetYaxis().UnZoom()
rectangle.Draw("same")
cnv.SaveAs(fnameout+"gaindrop_2d_anomaly_blackened.pdf")
cnv.SaveAs(fnameout+"gaindrop.pdf")

cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
h2.Draw("colz")
h2.GetYaxis().UnZoom()
ROOT.gPad.UnZoomed()
if(logz): ROOT.gPad.SetLogz()
h2Full.Draw("scat same")
h2Full.GetYaxis().UnZoom()
rectangle.Draw("same")
cnv.SaveAs(fnameout+"gaindrop_2d_fullarea.pdf")
cnv.SaveAs(fnameout+"gaindrop.pdf)")

print("=================== summary ===================")
print("avg in problematic blob=%g, rect=%g, rest=%g --> drop(blob)=%g, drop(rect)=%g for blob area fraction=%.1f%%" % (avg_blob, avg_rect, avg_rest, drop_blob, drop_rect, area_blob/area_full*100))
print("current of singular point=%g, rest=%g --> drop(singular)=%g%%" % (deepest, avg_rest, drop_singular))
