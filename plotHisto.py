#!/bin/python
from ROOT import gROOT, TCanvas, TF1, TFile
import ROOT

class HistogramFile(object):
	def __init__(self, filename):
		self.filename = filename

	def __enter__(self):
		self.file = ROOT.TFile.Open(self.filename, 'read')
		return self

	def __exit__(self, exception_type, exception_value, traceback):
		self.file.Close()

	def get_histogram(self, name):
		"""Return the histogram identified by name from the file.
		"""
		# TFile::Get() returns pointer to an object stored in a ROOT file.
		hist = self.file.Get(name)
		if hist:
			return hist
		else:
			raise RuntimeError('Unable to retrieve histogram named {0} from {1}'.format(name, self.filename))



#
f1 = TFile.Open('Hadro.root', 'read')
c1 = TCanvas('c1', 'Test Energy deposition')
c1.SetGridx()
c1.SetGridy()
print(f1.ls())
hist = f1.Get('h6')
hist.Draw()
c1.Update()
c1.SaveAs('plot.pdf')
f1.Close()
