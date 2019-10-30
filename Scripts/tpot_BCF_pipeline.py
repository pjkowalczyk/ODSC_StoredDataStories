import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from tpot.builtins import StackingEstimator
from sklearn.preprocessing import FunctionTransformer
from copy import copy

# NOTE: Make sure that the class is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR',
                        dtype=np.float64)
features = tpot_data.drop('target', axis=1).values
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'].values,
                             random_state=42)

# Score on the training set was:-0.2862459556813895
exported_pipeline = make_pipeline(
    make_union(
        FunctionTransformer(copy),
        StackingEstimator(estimator=ExtraTreesRegressor(bootstrap=True,
                                                        max_features=0.7000001,
                                                        min_samples_leaf=14,
                                                        min_samples_split=13,
                                                        n_estimators=100))
    ),
    ExtraTreesRegressor(bootstrap=False, max_features=0.8500000000000001,
                        min_samples_leaf=2, min_samples_split=7,
                        n_estimators=100)
)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
