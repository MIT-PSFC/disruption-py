#!/usr/bin/env python3

# TODO: https://github.com/MIT-PSFC/disruption-py/pull/271
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring

import multiprocessing
import threading
from enum import Enum
from typing import Callable

from disruption_py.config import config
from disruption_py.core.retrieval_manager import RetrievalManager
from disruption_py.inout.sql import ShotDatabase
from disruption_py.settings import OutputSetting, OutputSettingParams


class MarkCompleteEnum(Enum):
    """
    sentinel value for signifying that task queue is complete.
    """

    MARK_COMPLETE = "MARK_COMPLETE"


MARK_COMPLETE = MarkCompleteEnum.MARK_COMPLETE


class Consumer(multiprocessing.Process):
    def __init__(
        self,
        task_queue,
        result_queue,
        retrieval_manager_initializer: Callable[..., RetrievalManager],
    ):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.retrieval_manager_initializer = retrieval_manager_initializer

    def run(self):
        retrieval_manager: RetrievalManager = self.retrieval_manager_initializer()
        while True:
            next_task = self.task_queue.get()
            if next_task is MARK_COMPLETE:
                # Signal that the task is done and exit the loop
                self.task_queue.task_done()
                break

            shot_id, answer = next_task(retrieval_manager)

            self.task_queue.task_done()
            self.result_queue.put((shot_id, answer))


class ShotTask:
    def __init__(self, shot_id, retrieval_settings):
        self.shot_id = shot_id
        self.retrieval_settings = retrieval_settings

    def __call__(self, retrieval_manager: RetrievalManager):
        result = retrieval_manager.get_shot_data(
            shot_id=self.shot_id,
            retrieval_settings=self.retrieval_settings,
        )
        return self.shot_id, result

    def __str__(self):
        return f"Task for {self.shot_id}"


class MultiprocessingShotRetriever:
    """
    A class to run shot retrievals in parallel.
    """

    def __init__(
        self,
        database: ShotDatabase,
        output_setting: OutputSetting,
        retrieval_manager_initializer: Callable[..., RetrievalManager],
        tokamak,
        logger,
        num_processes=1,
    ):

        self.task_queue = multiprocessing.JoinableQueue()
        self.result_queue = multiprocessing.Queue()

        self.output_setting = output_setting
        self.database = database
        self.tokamak = tokamak
        self.logger = logger
        self.result_thread = threading.Thread(target=self._result_processor)

        self.consumers = [
            Consumer(
                task_queue=self.task_queue,
                result_queue=self.result_queue,
                retrieval_manager_initializer=retrieval_manager_initializer,
            )
            for _ in range(min(num_processes, config().MAX_PROCESSES))
        ]

    def _result_processor(self):
        while True:
            shot_id, result = self.result_queue.get()
            self.logger.info("Processing result for shot: %s", shot_id)
            if result is MARK_COMPLETE:
                break
            if result is None:
                self.logger.warning(
                    "Not outputting data for shot %s, data is None.", shot_id
                )
                continue
            self.output_setting.output_shot(
                OutputSettingParams(
                    shot_id=shot_id,
                    result=result,
                    database=self.database,
                    tokamak=self.tokamak,
                    logger=self.logger,
                )
            )

    def run(self, shotlist_list, retrieval_settings, await_complete=True):

        if not self.result_thread.is_alive():
            self.result_thread.start()

        for w in self.consumers:
            if not w.is_alive():
                w.start()

        for shot_id in shotlist_list:
            task = ShotTask(
                shot_id=shot_id,
                retrieval_settings=retrieval_settings,
            )
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

        # Signal the result processing thread to stop once completed processing
        # and wait for it to finish
        self.result_queue.put((None, MARK_COMPLETE))
        self.result_thread.join()
