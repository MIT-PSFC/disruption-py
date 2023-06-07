import pickle
import joblib
import numpy as np
import h5py 
from sklearn.ensemble import RandomForestClassifier


class HDF5RFModel():
    def __init__(self, n_trees, n_nodes, n_classes, value_arr, tree_start_arr, feature_arr, children_left_arr, children_right_arr, threshold_arr, stride):
        self.n_trees = n_trees
        self.n_nodes = n_nodes
        self.n_classes = n_classes
        self.value_arr = value_arr
        self.tree_start_arr = tree_start_arr
        self.feature_arr = feature_arr
        self.children_left_arr = children_left_arr
        self.children_right_arr = children_right_arr
        self.threshold_arr = threshold_arr
        self.stride = stride
        self.classes_ = [0, 1]

    def predict_proba(self, X):
        predictions = np.zeros((len(X)//self.stride,))
        results = np.zeros((self.n_trees,))
        for j in np.arange(len(X)//self.stride):
            X_eval = X[self.stride*j, :]
            for i in np.arange(self.n_trees):
                current_node = 0
                current_index = self.tree_start_arr[i]+current_node
                while self.children_left_arr[current_index] != self.children_right_arr[current_index]:
                    if X_eval[self.feature_arr[current_index]] <= self.threshold_arr[current_index]:
                        current_node = self.children_left_arr[current_index]
                    else:
                        current_node = self.children_right_arr[current_index]
                    current_index = self.tree_start_arr[i]+current_node
                results[i] = self.value_arr[current_index][1] / \
                    np.sum(self.value_arr[current_index])
            predictions[j] = np.sum(results)/self.n_trees
        return predictions


def load_model_from_hdf5(model_path, model_type='RandomForestClassifier'):
    if model_type == 'RandomForestClassifier':
        f_dict = h5py.File(model_path, 'r')
        n_trees = f_dict['n_trees'][0]
        n_nodes = f_dict['n_nodes'][0]
        n_classes = f_dict['n_classes'][0]
        value_arr = f_dict['value'][:]
        tree_start_arr = f_dict['tree_start'][:]
        feature_arr = f_dict['feature'][:]
        children_left_arr = f_dict['children_left'][:]
        children_right_arr = f_dict['children_right'][:]
        threshold_arr = f_dict['threshold'][:]
        stride = 1
        return HDF5RFModel(n_trees, n_nodes, n_classes, value_arr, tree_start_arr,
                           feature_arr, children_left_arr, children_right_arr, threshold_arr, stride)
    else:
        raise ValueError('Model type not supported')

def load_model(model_path):
    if model_path.endswith('.joblib'):
        model = joblib.load(model_path)
    elif model_path.endswith('.pkl'):
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
    elif model_path.endswith('.h5'):
        model = load_model_from_hdf5(model_path)
    else:
        raise ValueError(f"Unknown model file type: {model_path}")
    return model
  
def create_model(model_type):
    if model_type == 'random_forest':
        model = RandomForestClassifier(
            n_estimators=245, max_depth=15, random_state=0)
    else:
        raise NotImplementedError('Only random forest implemented')
    return model

