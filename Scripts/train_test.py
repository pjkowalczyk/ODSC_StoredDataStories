### Instantiate environment
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import numpy as np

#### AOH

df = pd.read_csv('cache/AOH_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogOH"])
y = df[["LogOH"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/AOH_Xtrain.csv", index = False)
X_test.to_csv("cache/AOH_Xtest.csv", index = False)
y_train.to_csv("cache/AOH_ytrain.csv", index = False)
y_test.to_csv("cache/AOH_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/AOH_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/AOH_Xtest_std.csv", index = False)

#### BCF

df = pd.read_csv('cache/BCF_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogBCF"])
y = df[["LogBCF"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/BCF_Xtrain.csv", index = False)
X_test.to_csv("cache/BCF_Xtest.csv", index = False)
y_train.to_csv("cache/BCF_ytrain.csv", index = False)
y_test.to_csv("cache/BCF_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/BCF_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/BCF_Xtest_std.csv", index = False)

#### BioHL

df = pd.read_csv('cache/BioHL_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogHalfLife"])
y = df[["LogHalfLife"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/BioHL_Xtrain.csv", index = False)
X_test.to_csv("cache/BioHL_Xtest.csv", index = False)
y_train.to_csv("cache/BioHL_ytrain.csv", index = False)
y_test.to_csv("cache/BioHL_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/BioHL_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/BioHL_Xtest_std.csv", index = False)

#### BP

df = pd.read_csv('cache/BP_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "BP"])
y = df[["BP"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/BP_Xtrain.csv", index = False)
X_test.to_csv("cache/BP_Xtest.csv", index = False)
y_train.to_csv("cache/BP_ytrain.csv", index = False)
y_test.to_csv("cache/BP_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/BP_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/BP_Xtest_std.csv", index = False)

#### HL

df = pd.read_csv('cache/HL_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogHL"])
y = df[["LogHL"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/HL_Xtrain.csv", index = False)
X_test.to_csv("cache/HL_Xtest.csv", index = False)
y_train.to_csv("cache/HL_ytrain.csv", index = False)
y_test.to_csv("cache/HL_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/HL_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/HL_Xtest_std.csv", index = False)

#### KM

df = pd.read_csv('cache/KM_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogKmHL"])
y = df[["LogKmHL"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/KM_Xtrain.csv", index = False)
X_test.to_csv("cache/KM_Xtest.csv", index = False)
y_train.to_csv("cache/KM_ytrain.csv", index = False)
y_test.to_csv("cache/KM_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/KM_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/KM_Xtest_std.csv", index = False)

#### KOA

df = pd.read_csv('cache/KOA_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogKOA"])
y = df[["LogKOA"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/KOA_Xtrain.csv", index = False)
X_test.to_csv("cache/KOA_Xtest.csv", index = False)
y_train.to_csv("cache/KOA_ytrain.csv", index = False)
y_test.to_csv("cache/KOA_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/KOA_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/KOA_Xtest_std.csv", index = False)

#### KOC

df = pd.read_csv('cache/KOC_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogKOC"])
y = df[["LogKOC"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/KOC_Xtrain.csv", index = False)
X_test.to_csv("cache/KOC_Xtest.csv", index = False)
y_train.to_csv("cache/KOC_ytrain.csv", index = False)
y_test.to_csv("cache/KOC_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/KOC_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/KOC_Xtest_std.csv", index = False)

#### logP

df = pd.read_csv('cache/logP_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogP"])
y = df[["LogP"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/logP_Xtrain.csv", index = False)
X_test.to_csv("cache/logP_Xtest.csv", index = False)
y_train.to_csv("cache/logP_ytrain.csv", index = False)
y_test.to_csv("cache/logP_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/logP_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/logP_Xtest_std.csv", index = False)

#### MP

df = pd.read_csv('cache/MP_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "MP"])
y = df[["MP"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/MP_Xtrain.csv", index = False)
X_test.to_csv("cache/MP_Xtest.csv", index = False)
y_train.to_csv("cache/MP_ytrain.csv", index = False)
y_test.to_csv("cache/MP_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/MP_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/MP_Xtest_std.csv", index = False)

#### RB

df = pd.read_csv('cache/RB_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "Ready_Biodeg"])
y = df[["Ready_Biodeg"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/RB_Xtrain.csv", index = False)
X_test.to_csv("cache/RB_Xtest.csv", index = False)
y_train.to_csv("cache/RB_ytrain.csv", index = False)
y_test.to_csv("cache/RB_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/RB_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/RB_Xtest_std.csv", index = False)

#### VP

df = pd.read_csv('cache/VP_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogVP"])
y = df[["LogVP"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/VP_Xtrain.csv", index = False)
X_test.to_csv("cache/VP_Xtest.csv", index = False)
y_train.to_csv("cache/VP_ytrain.csv", index = False)
y_test.to_csv("cache/VP_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/VP_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/VP_Xtest_std.csv", index = False)

#### WS

df = pd.read_csv('cache/WS_features.csv')
# df.sample(5).head()

### Training & Test Datasets
X = df.drop(columns=["Canonical_QSARr", "InChI_Code_QSARr", "LogMolar"])
y = df[["LogMolar"]]
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state = 42,
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

X_train.to_csv("cache/WS_Xtrain.csv", index = False)
X_test.to_csv("cache/WS_Xtest.csv", index = False)
y_train.to_csv("cache/WS_ytrain.csv", index = False)
y_test.to_csv("cache/WS_ytest.csv", index = False)
pd.DataFrame(X_train_standard).to_csv("cache/WS_Xtrain_std.csv", index = False)
pd.DataFrame(X_test_standard).to_csv("cache/WS_Xtest_std.csv", index = False)
