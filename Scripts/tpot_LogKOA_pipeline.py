import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.svm import LinearSVR
from tpot.builtins import StackingEstimator

# NOTE: Make sure that the class is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1).values
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'].values, random_state=42)

# Score on the training set was:-0.31986686085027194
exported_pipeline = make_pipeline(
    StackingEstimator(estimator=LinearSVR(C=0.01, dual=True, epsilon=0.0001, loss="epsilon_insensitive", tol=1e-05)),
    GradientBoostingRegressor(alpha=0.75, learning_rate=0.1, loss="huber", max_depth=9, max_features=0.6500000000000001, min_samples_leaf=2, min_samples_split=19, n_estimators=100, subsample=0.7000000000000001)
)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
