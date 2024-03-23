from enum import Enum
from typing import Callable, List, Dict
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
    def __init__(self, task_queue, result_queue, process_prop_initializers : Dict[str, Callable]):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.initialized_process_props = {prop : initializer() for prop, initializer in process_prop_initializers.items()}

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is MARK_COMPLETE:
                # Signal that the task is done and exit the loop
                self.task_queue.task_done()
                break

            shot_id, answer = next_task(self.initialized_process_props)
            
            self.task_queue.task_done()
            self.result_queue.put((shot_id, answer))
        return


class ShotTask:
    def __init__(self, shot_creator_f, shot_id, shot_settings):
        self.shot_creator_f = shot_creator_f
        self.shot_id = shot_id
        self.shot_settings = shot_settings

    def __call__(self, initialized_process_props: Dict):
        result = self.shot_creator_f(
            shot_id=self.shot_id, 
            shot_settings=self.shot_settings, 
            **initialized_process_props
        )
        return self.shot_id, result

    def __str__(self):
        return f'Task for {self.shot_id}'


class MultiprocessingShotRetriever:
    '''
    A class to run shot retrievals in parallel.
    '''
    
    def __init__(
        self, 
        shot_settings: ShotSettings, 
        database : ShotDatabase, 
        output_type_request: OutputTypeRequest, 
        process_prop_initializers : Dict[str, Callable],
        tokamak, 
        logger, 
        num_processes=1, 
    ):
        
        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.shot_settings = shot_settings
        self.output_type_request = output_type_request
        self.database = database
        self.tokamak = tokamak
        self.logger = logger
        self.result_thread = threading.Thread(target=self._result_processor)

        self.consumers = [
            Consumer(self.task_queue, self.result_queue, process_prop_initializers=process_prop_initializers) 
            for _ in range(min(num_processes, MAX_PROCESSES))
        ]
  
    def _result_processor(self):
        while True:
            shot_id, result = self.result_queue.get()
            self.logger.info(f"Processing result for shot: {shot_id}")
            if result is MARK_COMPLETE:
                break
            elif result is None:
                self.logger.warning(f"Not outputting data for shot {shot_id}, data is None.")
                continue
            self.output_type_request.output_shot(
                ResultOutputTypeRequestParams(
                    shot_id=shot_id,
                    result=result, 
                    database=self.database, 
                    tokamak=self.tokamak, 
                    logger=self.logger,
                )
            )

    def run(self, shot_creator_f, shot_ids_list, await_complete=True):
        
        if not self.result_thread.is_alive():
            self.result_thread.start()
   
        for w in self.consumers:
            if not w.is_alive():
                w.start()
        
        for shot_id in shot_ids_list:
            task = ShotTask(shot_creator_f=shot_creator_f, shot_id=shot_id, shot_settings=self.shot_settings)
            self.task_queue.put(task)
   
        if await_complete:
            self.await_complete()

        return True

    def await_complete(self):
         # Signal the consumers to stop once completed processing and wait for it to finish
        for _ in self.consumers:
            self.task_queue.put(MARK_COMPLETE)
        self.task_queue.join()
        
        for consumer in self.consumers:
            consumer.join()

        # Signal the result processing thread to stop once completed processing and wait for it to finish
        self.result_queue.put((None, MARK_COMPLETE))
        self.result_thread.join()