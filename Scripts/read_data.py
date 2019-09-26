### Instantiate environment
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd

### Read data

#### AOH: atmospheric hydroxylation rate
train_df = PandasTools.LoadSDF("data/TR_AOH_516.sdf")
test_df = PandasTools.LoadSDF("data/TST_AOH_176.sdf")
AOH = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogOH"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogOH"]]], ignore_index = True)
AOH['LogOH'] = pd.to_numeric(AOH['LogOH'])

for index, row in AOH.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        AOH.loc[index, 'molecule'] = True
    else:
        AOH.loc[index, 'molecule'] = False
AOH = AOH[AOH['molecule'] == True]
AOH = AOH.drop(columns="molecule")

AOH.to_csv("cache/AOH.csv", index = False)

#### BCF: bioconcentration factor
train_df = PandasTools.LoadSDF("data/TR_BCF_469.sdf")
test_df = PandasTools.LoadSDF("data/TST_BCF_157.sdf")
BCF = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogBCF"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogBCF"]]], ignore_index = True)
BCF['LogBCF'] = pd.to_numeric(BCF['LogBCF'])

for index, row in BCF.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        BCF.loc[index, 'molecule'] = True
    else:
        BCF.loc[index, 'molecule'] = False
BCF = BCF[BCF['molecule'] == True]
BCF = BCF.drop(columns="molecule")

BCF.to_csv("cache/BCF.csv", index = False)

#### BioHL: biodegradability half-life
train_df = PandasTools.LoadSDF("data/TR_BioHL_112.sdf")
test_df = PandasTools.LoadSDF("data/TST_BioHL_38.sdf")
BioHL = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogHalfLife"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogHalfLife"]]], ignore_index = True)
BioHL['LogHalfLife'] = pd.to_numeric(BioHL['LogHalfLife'])

for index, row in BioHL.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        BioHL.loc[index, 'molecule'] = True
    else:
        BioHL.loc[index, 'molecule'] = False
BioHL = BioHL[BioHL['molecule'] == True]
BioHL = BioHL.drop(columns="molecule")

BioHL.to_csv("cache/BioHL.csv", index = False)

#### BP: boiling point
train_df = PandasTools.LoadSDF("data/TR_BP_4077.sdf")
test_df = PandasTools.LoadSDF("data/TST_BP_1358.sdf")
BP = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "BP"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "BP"]]], ignore_index = True)
BP['BP'] = pd.to_numeric(BP['BP'])

for index, row in BP.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        BP.loc[index, 'molecule'] = True
    else:
        BP.loc[index, 'molecule'] = False
BP = BP[BP['molecule'] == True]
BP = BP.drop(columns="molecule")

BP.to_csv("cache/BP.csv", index = False)

#### HL: Henry's Law constant
train_df = PandasTools.LoadSDF("data/TR_HL_441.sdf")
test_df = PandasTools.LoadSDF("data/TST_HL_150.sdf")
HL = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogHL"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogHL"]]], ignore_index = True)
HL['LogHL'] = pd.to_numeric(HL['LogHL'])

for index, row in HL.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        HL.loc[index, 'molecule'] = True
    else:
        HL.loc[index, 'molecule'] = False
HL = HL[HL['molecule'] == True]
HL = HL.drop(columns="molecule")

HL.to_csv("cache/HL.csv", index = False)

#### KM: fish biotransformation half-life
train_df = PandasTools.LoadSDF("data/TR_KM_405.sdf")
test_df = PandasTools.LoadSDF("data/TST_KM_136.sdf")
KM = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKmHL"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKmHL"]]], ignore_index = True)
KM['LogKmHL'] = pd.to_numeric(KM['LogKmHL'])

for index, row in KM.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        KM.loc[index, 'molecule'] = True
    else:
        KM.loc[index, 'molecule'] = False
KM = KM[KM['molecule'] == True]
KM = KM.drop(columns="molecule")

KM.to_csv("cache/KM.csv", index = False)

#### KOA: octanol-air partition coefficient
train_df = PandasTools.LoadSDF("data/TR_KOA_202.sdf")
test_df = PandasTools.LoadSDF("data/TST_KOA_68.sdf")
KOA = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKOA"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKOA"]]], ignore_index = True)
KOA['LogKOA'] = pd.to_numeric(KOA['LogKOA'])

for index, row in KOA.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        KOA.loc[index, 'molecule'] = True
    else:
        KOA.loc[index, 'molecule'] = False
KOA = KOA[KOA['molecule'] == True]
KOA = KOA.drop(columns="molecule")

KOA.to_csv("cache/KOA.csv", index = False)

#### KOC: soil adsorption coefficient
train_df = PandasTools.LoadSDF("data/TR_KOC_545.sdf")
test_df = PandasTools.LoadSDF("data/TST_KOC_184.sdf")
KOC = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKOC"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogKOC"]]], ignore_index = True)
KOC['LogKOC'] = pd.to_numeric(KOC['LogKOC'])

for index, row in KOC.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        KOC.loc[index, 'molecule'] = True
    else:
        KOC.loc[index, 'molecule'] = False
KOC = KOC[KOC['molecule'] == True]
KOC = KOC.drop(columns="molecule")

KOC.to_csv("cache/KOC.csv", index = False)

#### logP: octanol-water partition coefficient
train_df = PandasTools.LoadSDF("data/TR_LogP_10537.sdf")
test_df = PandasTools.LoadSDF("data/TST_LogP_3513.sdf")
logP = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogP"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogP"]]], ignore_index = True)
logP['LogP'] = pd.to_numeric(logP['LogP'])

for index, row in logP.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        logP.loc[index, 'molecule'] = True
    else:
        logP.loc[index, 'molecule'] = False
logP = logP[logP['molecule'] == True]
logP = logP.drop(columns="molecule")

logP.to_csv("cache/logP.csv", index = False)

#### MP: melting point
train_df = PandasTools.LoadSDF("data/TR_MP_6486.sdf")
test_df = PandasTools.LoadSDF("data/TST_MP_2167.sdf")
MP = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "MP"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "MP"]]], ignore_index = True)
MP['MP'] = pd.to_numeric(MP['MP'])

for index, row in MP.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        MP.loc[index, 'molecule'] = True
    else:
        MP.loc[index, 'molecule'] = False
MP = MP[MP['molecule'] == True]
MP = MP.drop(columns="molecule")

MP.to_csv("cache/MP.csv", index = False)

#### RB: readily biodegradable
train_df = PandasTools.LoadSDF("data/TR_RBioDeg_1197.sdf")
test_df = PandasTools.LoadSDF("data/TST_RBioDeg_411.sdf")
RB = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "Ready_Biodeg"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "Ready_Biodeg"]]], ignore_index = True)
RB['Ready_Biodeg'] = pd.to_numeric(RB['Ready_Biodeg'])

for index, row in RB.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        RB.loc[index, 'molecule'] = True
    else:
        RB.loc[index, 'molecule'] = False
RB = RB[RB['molecule'] == True]
RB = RB.drop(columns="molecule")

RB.to_csv("cache/RB.csv", index = False)

#### VP: vapor pressure
train_df = PandasTools.LoadSDF("data/TR_VP_2034.sdf")
test_df = PandasTools.LoadSDF("data/TST_VP_679.sdf")
VP = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogVP"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogVP"]]], ignore_index = True)
VP['LogKOC'] = pd.to_numeric(VP['LogVP'])

for index, row in VP.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        VP.loc[index, 'molecule'] = True
    else:
        VP.loc[index, 'molecule'] = False
VP = VP[VP['molecule'] == True]
VP = VP.drop(columns="molecule")

VP.to_csv("cache/VP.csv", index = False)

#### WS: water solubility
train_df = PandasTools.LoadSDF("data/TR_WS_3158.sdf")
test_df = PandasTools.LoadSDF("data/TST_WS_1066.sdf")
WS = pd.concat([train_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogMolar"]],
                 test_df[["Canonical_QSARr", "InChI_Code_QSARr", "LogMolar"]]], ignore_index = True)
WS['LogMolar'] = pd.to_numeric(WS['LogMolar'])

for index, row in WS.iterrows():
    m = Chem.MolFromSmiles(row['Canonical_QSARr'])
    if m is not None:
        WS.loc[index, 'molecule'] = True
    else:
       WS.loc[index, 'molecule'] = False
WS = WS[WS['molecule'] == True]
WS = WS.drop(columns="molecule")

WS.to_csv("cache/WS.csv", index = False)

