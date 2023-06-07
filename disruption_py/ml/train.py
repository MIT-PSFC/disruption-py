from sklearn.ensemble import RandomForestClassifier
import numpy as np

def grid_search(x_train, y_train, x_test, y_test, model_type, grid, **kwargs):
    # Grid is a dictionary where each key is a parameter and each value is a list of values to try
    # kwargs are any other parameters to pass to the model
    searches = []
    for param in grid:
        for val in grid[param]:
            if len(searches) == 0:
                searches.append({param: val})
            else:
                searches = [dict(search, **{param: val})
                            for search in searches]
    results = []
    for search in searches:
        result = train_local(x_train, y_train, x_test,
                             y_test,  **search, **kwargs)
        results.append({'params': search, 'result': result})
    return results

def train_local(x_train, y_train, x_test, y_test, **kwargs):
    if 'omit_after_quench' in kwargs and kwargs['omit_after_quench'] is not None:
        raise NotImplementedError("Talk to Cristina about implementing")
    if 'features' in kwargs and kwargs['features'] is not None:
        features = kwargs['features']
        x_train = x_train[features]
    if 'model' in kwargs and kwargs['model'] is not None:
        raise NotImplementedError('Only random forest implemented')
    else:
        n_estimators = kwargs.pop('n_estimators', 245)
        max_depth = kwargs.pop('max_depth', 15)
        random_state = kwargs.pop('random_state', 0)
        n_jobs = kwargs.pop('n_jobs', -1)
        model = RandomForestClassifier(
            n_estimators=n_estimators, max_depth=max_depth, random_state=random_state, n_jobs=n_jobs, **kwargs)
        model.fit(x_train, y_train)
        return {'model': model,
                'predictions': model.predict(x_train),
                'predictions_proba': model.predict_proba(x_train),
                'features': model.feature_importances_,
                'std_imp': np.std([tree.feature_importances_ for tree in model.estimators_], axis=0),
                'score': model.score(x_test, y_test)}
