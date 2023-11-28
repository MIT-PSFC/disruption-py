from typing import List, Dict
from disruption_py.settings.output_type_requests import OutputTypeRequestParams, output_type_request_runner
from disruption_py.utils.mappings.tokemak import Tokemak
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
    def __init__(self, shot_creator_f, shot_args: Dict):
        self.shot_creator_f = shot_creator_f
        self.shot_args = shot_args

    def __call__(self, task_args: Dict):
        result = self.shot_creator_f(**task_args, **self.shot_args)
        return result

    def __str__(self):
        return f'Task for {self.shot_args.get("shot_id", "Unknown shot number")}'


class MultiprocessingShotRetriever:
    '''
    A class to run shot retrievals in parallel.
    '''
    
    def __init__(self, output_type_request, tokemak, logger, num_processes=8, database_initializer_f = None):
        
        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.output_type_request = output_type_request
        self.tokemak = tokemak
        self.logger = logger
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
            output_type_request_runner(self.output_type_request, OutputTypeRequestParams(result, self.tokemak, self.logger))

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

        self.output_type_request.stream_output_cleanup()
        return self.output_type_request.get_results()
