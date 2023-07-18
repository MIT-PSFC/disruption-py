from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import make_pipeline
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
    omit_after_quench = kwargs.pop('omit_after_quench', None)
    features = kwargs.pop('features', None)
    model = kwargs.pop('model', "RandomForestClassifier")
    class_weight = kwargs.pop('class_weight', 'balanced')
    if omit_after_quench is not None:
        raise NotImplementedError("Talk to Cristina about implementing")
    if features is not None:
        x_train = x_train[features]

    if model == "LogisticRegression":
        C = kwargs.pop('C', 1.0)
        random_state = kwargs.pop('random_state', 0)
        log_reg= LogisticRegression(C=C, random_state=random_state, max_iter=10000,class_weight='balanced', **kwargs)
        poly = PolynomialFeatures(degree=4)
        scaler = StandardScaler()
        model = make_pipeline(scaler, poly, log_reg)
        model.fit(x_train, y_train)
        return {'model': model,
                'predictions': model.predict(x_train),
                'predictions_proba': model.predict_proba(x_train),
                'features': model.named_steps['logisticregression'].coef_[0],  # If binary classification
                'score': model.score(x_test, y_test)}
    elif model == "SupportVectorMachine":
        C = kwargs.pop('C', 1.0)
        random_state = kwargs.pop('random_state', 0)
        model = SVC(C=C, random_state=random_state, kernel='rbf', probability=True, class_weight='balanced', **kwargs)
        model.fit(x_train, y_train)
        feature_importances = [0]*len(x_train.columns)    
        return {'model': model,
                'predictions': model.predict(x_train),
                'predictions_proba': model.predict_proba(x_train),
                'features': feature_importances,   # Placeholder for SVM, no direct feature importance
                'score': model.score(x_test, y_test)}
    elif model == "RandomForestClassifier":
        n_estimators = kwargs.pop('n_estimators', 245)
        max_depth = kwargs.pop('max_depth', 15)
        random_state = kwargs.pop('random_state', 0)
        n_jobs = kwargs.pop('n_jobs', -1)
        model = RandomForestClassifier(
            n_estimators=n_estimators, max_depth=max_depth, random_state=random_state, n_jobs=n_jobs,class_weight=class_weight,**kwargs)
        model.fit(x_train, y_train)
        return {'model': model,
                'predictions': model.predict(x_train),
                'predictions_proba': model.predict_proba(x_train),
                'features': model.feature_importances_,
                'std_imp': np.std([tree.feature_importances_ for tree in model.estimators_], axis=0),
                'score': model.score(x_test, y_test)}  
    else:
        raise NotImplementedError('Only random forest, logistic regression, and support vector machines are implemented')
