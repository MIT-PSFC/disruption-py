from enum import Enum
from typing import List, Dict
from disruption_py.databases.database import ShotDatabase
from disruption_py.settings import ShotSettings, OutputTypeRequest, ResultOutputTypeRequestParams, FinishOutputTypeRequestParams
from disruption_py.utils.constants import MAX_PROCESSES
import multiprocessing
import threading

# define a sentinel value for signifying that task queue is complete
class MarkCompleteEnum(Enum):
    MarkComplete = 'MarkComplete'
MARK_COMPLETE = MarkCompleteEnum.MarkComplete

class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, initialize_database=True, database_initializer_f = None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self._consumer_database = None

        if initialize_database and database_initializer_f is None:
            raise ValueError("Must provide database_initializer_f if initialize_database is True")
        self.initialize_database = initialize_database
        self.database_initializer_f = database_initializer_f

    @property
    def consumer_database(self):
        if self._consumer_database is None:
            self._consumer_database = self.database_initializer_f()
        return self._consumer_database

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is MARK_COMPLETE:
                # Signal that the task is done and exit the loop
                self.task_queue.task_done()
                break

            task_args = {}
            if self.initialize_database:
                task_args['sql_database'] = self.consumer_database
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
    
    def __init__(self, shot_settings: ShotSettings, database : ShotDatabase, output_type_request: OutputTypeRequest, tokamak, logger, num_processes=8, database_initializer_f = None):
        
        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.shot_settings = shot_settings
        self.output_type_request = output_type_request
        self.database = database
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
            if result is MARK_COMPLETE:
                break
            elif result is None:
                continue
            self.output_type_request.output_shot(ResultOutputTypeRequestParams(result, self.database, self.tokamak, self.logger))

    def run(self, shot_creator_f, shot_ids_list, should_finish=True):
        
        if not self.result_thread.is_alive():
            self.result_thread.start()
   
        for w in self.consumers:
            if not w.is_alive():
                w.start()
        
        for shot_id in shot_ids_list:
            task = ShotTask(shot_creator_f=shot_creator_f, shot_id=shot_id, shot_settings=self.shot_settings)
            self.task_queue.put(task)
   
        if should_finish:
            return self.finish_and_await_results()

        return True

    def finish_and_await_results(self):
         # Signal the consumers to stop once completed processing and wait for it to finish
        for _ in self.consumers:
            self.task_queue.put(MARK_COMPLETE)
        self.task_queue.join()
        
        # Signal the result processing thread to stop once completed processing and wait for it to finish
        self.result_queue.put(MARK_COMPLETE)
        self.result_thread.join()

        finish_output_type_request_params = FinishOutputTypeRequestParams(self.tokamak, self.logger)
        result = self.output_type_request.get_results(finish_output_type_request_params)
        self.output_type_request.stream_output_cleanup(finish_output_type_request_params)
        return result
