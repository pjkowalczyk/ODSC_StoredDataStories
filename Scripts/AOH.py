ts# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

### Instantiate environment
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import PandasTools
import pandas as pd
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 
import numpy as np
import math
from sklearn.ensemble import RandomForestRegressor

### Read data
train_df = PandasTools.LoadSDF("data/TR_AOH_516.sdf")
test_df = PandasTools.LoadSDF("data/TST_AOH_176.sdf")

### Concatenate data
AOH = pd.concat([train_df[["Canonical_QSARr", "LogOH"]],
                 test_df[["Canonical_QSARr", "LogOH"]]], ignore_index = True)
AOH['LogOH'] = pd.to_numeric(AOH['LogOH'])

### Calculate features
nms = [x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
for i in range(len(AOH)):
    descrs = calc.CalcDescriptors(Chem.MolFromSmiles(AOH.iloc[i, 0]))
    for x in range(len(descrs)):
        AOH.at[i, str(nms[x])] = descrs[x]
AOH = AOH.dropna()

AOH.shape

### Training & Test Datasets
X = AOH.drop(columns=["Canonical_QSARr", "LogOH"])
y = AOH[["LogOH"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 350,
                                                    test_size = 0.2)

### Identify / remove near-zero variance descriptors
def variance_threshold_selector(data, threshold = 0.5):
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    return data[data.columns[selector.get_support(indices = True)]]

nzv = variance_threshold_selector(X_train, 0.0)

X_train = X_train[nzv.columns]
X_test = X_test[nzv.columns]

### Identify / remove highly correlated descriptors
corr_matrix = X_train.corr().abs()
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape),
                                  k = 1).astype(np.bool))
to_drop = [column for column in upper.columns
           if any(upper[column] > 0.85)]

X_train = X_train[X_train.columns.drop(to_drop)]
X_test = X_test[X_test.columns.drop(to_drop)]

### standardize features by removing the mean and scaling to unit variance
scaler = StandardScaler()
scaler.fit(X_train)

X_train_standard = scaler.transform(X_train)
X_test_standard = scaler.transform(X_test)

# ####
# #### TPOT
# ####

from tpot import TPOTRegressor
tpot = TPOTRegressor(generations=10, population_size=50, verbosity=2)
tpot.fit(X_train_standard, y_train)
print(tpot.score(X_test_standard, y_test))
tpot.export('tpot_AOH_pipeline.py')

# ####
# #### Use best pipeline
# ####

from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import LassoLarsCV
from sklearn.pipeline import make_pipeline, make_union
from sklearn.svm import LinearSVR
from tpot.builtins import StackingEstimator

### dataset for prediction modeling
dframe = PandasTools.LoadSDF("data/TST_AOH_176.sdf")

dframe = dframe[["Canonical_QSARr", "LogOH"]]

### Calculate features
nms = [x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
for i in range(len(dframe)):
    descrs = calc.CalcDescriptors(Chem.MolFromSmiles(dframe.iloc[i, 0]))
    for x in range(len(descrs)):
        dframe.at[i, str(nms[x])] = descrs[x]
dframe = dframe.dropna()

dframe.shape

dframe['LogOH'] = pd.to_numeric(dframe['LogOH'])
observed = np.array(dframe["LogOH"])
features = dframe[dframe.columns.drop("Canonical_QSARr", "LogOH")]
features = features[nzv.columns]
features = features[features.columns.drop(to_drop)]
features_standard = scaler.transform(features)

# NOTE: Make sure that the class is labeled 'target' in the data file
# tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR',
#                         dtype=np.float64)
# features = tpot_data.drop('target', axis=1).values
# training_features, testing_features, training_target, testing_target = \
#             train_test_split(features, tpot_data['target'].values,
#                              random_state=42)

training_features = X_train_standard
testing_features = X_test_standard
training_target = np.array(y_train['LogOH'])
testing_target =  np.array(y_test['LogOH'])

# Score on the training set was:-0.23600217676872473
exported_pipeline = make_pipeline(
    StackingEstimator(estimator=LassoLarsCV(normalize=True)),
    StackingEstimator(estimator=LinearSVR(C=10.0, dual=True, epsilon=0.0001,
                                          loss="epsilon_insensitive",
                                          tol=0.1)),
    ExtraTreesRegressor(bootstrap=False, max_features=0.9000000000000001,
                        min_samples_leaf=1, min_samples_split=9,
                        n_estimators=100)
)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)

#.....#####.....#####.....#

obs = np.array(y_test["LogOH"])
import matplotlib.pyplot as plt
plt.scatter(obs, results)
plt.ylabel('some numbers')
plt.show()

results1 = exported_pipeline.predict(training_features)
obs = np.array(y_train["LogOH"])
import matplotlib.pyplot as plt
plt.scatter(obs, results1)
plt.ylabel('some numbers')
plt.show()

results2 = exported_pipeline.predict(features_standard)
import matplotlib.pyplot as plt
plt.scatter(observed, results2)
plt.ylabel('some numbers')
plt.show()


