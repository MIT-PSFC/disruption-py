import multiprocessing
import threading
import pandas as pd

class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, initialize_database=True, database_initializer_f = None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self._sql_database = None

        if initialize_database and database_initializer_f is None:
            raise ValueError("Must provide database_initializer_f if initialize_database is True")
        self.initialize_database = initialize_database
        self.database_initializer_f = database_initializer_f

    @property
    def sql_database(self):
        if self._sql_database is None:
            self._sql_database = self.database_initializer_f()
        return self._sql_database

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Signal that the task is done and exit the loop
                self.task_queue.task_done()
                break

            task_args = {}
            if self.initialize_database:
                task_args['sql_database'] = self.sql_database
            answer = next_task(**task_args)
   
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return


class ShotTask:
    def __init__(self, shot_creator_f, shot_args):
        self.shot_creator_f =shot_creator_f
        self.shot_args = shot_args

    def __call__(self, **task_args):
        result = self.shot_creator_f(**task_args, **self.shot_args)
        return result

    def __str__(self):
        return f'Task for {self.shot_args.get("shot_id", "Unknown shot number")}'


# Functions used to stream the output of the shot retrievals
def get_hdf5_stream_ouput_fs(filepath):
    
    store = pd.HDFStore('streamed_data.h5', mode='w')
    def ouput_shot_df_hdf5(result, results):
        if results is None:
            results = 0
        store.append(f'df_{results}', result, format='table', data_columns=True)
        return results + 1
    
    def stream_output_cleanup_hdf5(results):
        store.close()
    
    return ouput_shot_df_hdf5, stream_output_cleanup_hdf5

def ouput_shot_df_list(result, results):
    if results is None:
        results = []
    results.append(result)
    return results

class MultiprocessingShotRetriever:
    '''
    A class to run shot retrievals in parallel.
    '''
    
    def __init__(self, num_processes=8, stream_ouput_process_f=ouput_shot_df_list, stream_ouput_cleanup_f=None, database_initializer_f = None):
        
        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.stream_ouput_process_f = stream_ouput_process_f
        self.stream_ouput_cleanup_f = stream_ouput_cleanup_f
        # Start the result processing thread
        self.results = None
        self.result_thread = threading.Thread(target=self._result_processor)

        self.consumers = [
            Consumer(self.task_queue, self.result_queue, database_initializer_f=database_initializer_f) 
            for _ in range(num_processes)
        ]
  
    def _result_processor(self):
        while True:
            result = self.result_queue.get()
            if result is None:
                break
            self.results = self.stream_ouput_process_f(result, self.results)

    def run(self, shot_creator_f, shot_id_list, shot_args_dict, should_finish=True):
        
        if not self.result_thread.is_alive():
            self.result_thread.start()
   
        for w in self.consumers:
            if not w.is_alive():
                w.start()
        
        for shot_id in shot_id_list:
            task = ShotTask(shot_creator_f, {**shot_args_dict, 'shot_id': shot_id})
            self.task_queue.put(task)
   
        if should_finish:
            return self.finish_and_await_results()

        return True

    def finish_and_await_results(self):
         # Signal the consumers to stop once completed processing and wait for it to finish
        for _ in self.consumers:
            self.task_queue.put(None)
        self.task_queue.join()
        
        # Signal the result processing thread to stop once completed processing and wait for it to finish
        self.result_queue.put(None)
        self.result_thread.join()

        if self.stream_ouput_cleanup_f is not None:
            self.stream_ouput_cleanup_f(self.results)
        
        return self.results
