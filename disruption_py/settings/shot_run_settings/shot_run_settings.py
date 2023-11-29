from dataclasses import dataclass, field
from typing import List
from disruption_py.settings.shot_run_settings.shot_data_requests import ShotDataRequest
from disruption_py.utils.mappings.mappings_helpers import map_string_attributes_to_enum

def default_tags():
    return ["all"]

@dataclass
class ShotRunSettings:
    run_methods : List[str] = field(default_factory=list)
    run_tags : List[str] = field(default_factory=default_tags)
    shot_data_requests : List[ShotDataRequest] = field(default_factory=list)
