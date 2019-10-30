import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor
from sklearn.linear_model import LassoLarsCV
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import StandardScaler
from tpot.builtins import StackingEstimator

# NOTE: Make sure that the class is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1).values
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'].values, random_state=42)

# Score on the training set was:-0.2956849979330415
exported_pipeline = make_pipeline(
    StandardScaler(),
    StackingEstimator(estimator=LassoLarsCV(normalize=False)),
    StackingEstimator(estimator=GradientBoostingRegressor(alpha=0.95, learning_rate=0.1, loss="ls", max_depth=5, max_features=0.9000000000000001, min_samples_leaf=20, min_samples_split=9, n_estimators=100, subsample=0.25)),
    ExtraTreesRegressor(bootstrap=False, max_features=0.3, min_samples_leaf=1, min_samples_split=5, n_estimators=100)
)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
