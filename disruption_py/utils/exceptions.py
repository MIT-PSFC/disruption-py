from disruption_py.machine.tokamak import Tokamak


class TokamakNotSupportedError(Exception):
    def __init__(self, tokamak: Tokamak, use_case: str):
        super().__init__(f"Tokamak {tokamak.value} not not supported for {use_case}.")
