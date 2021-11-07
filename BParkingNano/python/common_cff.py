import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

K_MASS = 0.493677
PI_MASS = 0.139571

def ufloat(expr, precision = -1, doc = ''):
  return Var('userFloat("%s")' % expr, 
             float, precision = precision, doc = doc)

def uint(expr, doc = ''):
  return Var('userInt("%s")' % expr, int, doc = doc)

def ubool(expr, precision = -1, doc = ''):
  return Var('userInt("%s") == 1' % expr, bool, doc = doc)

