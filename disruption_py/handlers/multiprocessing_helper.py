from typing import List, Dict
from disruption_py.settings import ShotSettings, ResultOutputTypeRequestParams, FinishOutputTypeRequestParams
from disruption_py.utils.constants import MAX_PROCESSES
from disruption_py.utils.mappings.tokamak import Tokamak
import multiprocessing
import threading

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
            answer = next_task(task_args)
   
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return


class ShotTask:
    def __init__(self, shot_creator_f, shot_id, shot_settings):
        self.shot_creator_f = shot_creator_f
        self.shot_id = shot_id
        self.shot_settings = shot_settings

    def __call__(self, task_args: Dict):
        result = self.shot_creator_f(shot_id=self.shot_id, shot_settings=self.shot_settings, **task_args)
        return result

    def __str__(self):
        return f'Task for {self.shot_id}'


class MultiprocessingShotRetriever:
    '''
    A class to run shot retrievals in parallel.
    '''
    
    def __init__(self, shot_settings: ShotSettings, tokamak, logger, num_processes=8, database_initializer_f = None):
        
        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.shot_settings = shot_settings
        self.tokamak = tokamak
        self.logger = logger
        self.result_thread = threading.Thread(target=self._result_processor)

        self.consumers = [
            Consumer(self.task_queue, self.result_queue, database_initializer_f=database_initializer_f) 
            for _ in range(min(num_processes, MAX_PROCESSES))
        ]
  
    def _result_processor(self):
        while True:
            result = self.result_queue.get()
            if result is None:
                break
            self.shot_settings.output_type_request.output_shot(ResultOutputTypeRequestParams(result, self.tokamak, self.logger))

    def run(self, shot_creator_f, shot_id_list, should_finish=True):
        
        if not self.result_thread.is_alive():
            self.result_thread.start()
   
        for w in self.consumers:
            if not w.is_alive():
                w.start()
        
        for shot_id in shot_id_list:
            task = ShotTask(shot_creator_f=shot_creator_f, shot_id=shot_id, shot_settings=self.shot_settings)
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

        finish_output_type_request_params = FinishOutputTypeRequestParams(self.tokamak, self.logger)
        self.shot_settings.output_type_request.stream_output_cleanup(finish_output_type_request_params)
        return self.shot_settings.output_type_request.get_results(finish_output_type_request_params)
