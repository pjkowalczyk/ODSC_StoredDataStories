# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 1 2018

@author: P J Kowalczyk
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import pandas as pd
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import MACCSkeys
from rdkit.RDLogger import logger
logger = logger()

# nms = names of calculated features
nms=[x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)

frame = PandasTools.LoadSDF('data/TST_WS_1066.sdf') # edit filename
w = Chem.SDWriter('cache/TST_WS_1066_descrs.sdf') # edit filename
for row in frame.itertuples():
    mol = Chem.MolFromSmiles(row.Canonical_QSARr)
    if mol is None: continue
    descrs = calc.CalcDescriptors(mol)
    for nm, v in zip(nms, descrs):
        mol.SetProp(nm, str(v))
        mol.SetProp('CAS', row.CAS)
        mol.SetProp('LogMolar', row.LogMolar) # edit target
        mol.SetProp('SMILES', row.Canonical_QSARr)
    w.write(mol)
 
frame2 = PandasTools.LoadSDF('cache/TST_WS_1066_descrs.sdf') # edit filename

frame2.to_csv('cache/TST_WS_1066_descrs.csv') # edit filename
