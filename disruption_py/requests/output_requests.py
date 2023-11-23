import pandas as pd
from abc import ABC, abstractmethod


class OutputProcessor(ABC):
    
    @abstractmethod
    def ouput_shot(self, result):
        pass
    
    def stream_output_cleanup(self):
        pass
    
    @abstractmethod
    def get_results(self):
        pass


class ListOuptutProcessor(OutputProcessor):
    def __init__(self):
        self.results = []
        
    def ouput_shot(self, result):
        self.results.append(result)
    
    def get_results(self):
        return self.results
    
    
class HDF5OutputProcessor(OutputProcessor):
    def __init__(self, filepath):
        self.filepath = filepath
        self.store = pd.HDFStore(filepath, mode='w')
        self.output_shot_count = 0

    def ouput_shot(self, result: pd.DataFrame):
        shot_id = result['shot'].iloc[0] if (not result.empty or result not in result.columns) else self.output_shot_count
        self.store.append(f'df_{shot_id}', result, format='table', data_columns=True)
        self.output_shot_count += 1
    
    def stream_output_cleanup(self):
        self.store.close()
    
    def get_results(self):
        return self.output_shot_count
