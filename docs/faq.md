## Installation { .doc .doc-heading }

### Issues with package versioning { .doc .doc-heading }
If the machine that you are working on has an outdated version of python, you may be unable to install.

#### On CMod { .doc .doc-heading }
On the mfe workstations use mferws02 or mferws03

### How can I use disruption_py if I don't have access to the GitHub Repository? { .doc .doc-heading }
If you are unable to access the GitHub repository you can manually install the package. Note that you may be installing an older version of DisruptionPy.

#### On CMod { .doc .doc-heading }
You can do this by running:
```bash
pip install /home/joshlor/disruption-py
```

### Stuck on Poetry Install { .doc .doc-heading }
Poetry might hang instead of asking to unlock the keyring, which is a [known bug](https://github.com/python-poetry/poetry/issues/8623).
As a temporary workaround while the bug is addressed, please set:
```bash
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
```

## Usage { .doc .doc-heading }

### How do I debug issues with DisruptionPy? { .doc .doc-heading }
To better understand the issue that you are seeing in DisruptionPy you can enable more detailed logging by passing `log_settings` to the `retrieval_settings` for your program. For more information on customizing logging, please see [`LogSettings`][disruption_py.settings.log_settings.LogSettings].

```python
import logging
from disruption_py.settings.retrieval_settings import ShotSettings
from disruption_py.settings.log_settings import LogSettings

ShotSettings(
    ...

    log_settings=LogSettings(
        log_to_console=False, # tell DisruptionPy to not log to the console
        log_file_path="path/to/log/file",
        file_log_level=logging.DEBUG, # tell DisruptionPy to log
    ),

    ...
)
```

### How can I get information on data accuracy? { .doc .doc-heading }
DisruptionPy provides the [`disruption_py run evaluate`][disruption_py-run-evaluate] command to provide insight on the current accuracy of DisruptionPy methods.

### Why does DisruptionPy log `%MDSPLUS-E-ERROR` for all shots after a certain shot number? { .doc .doc-heading }
You may have a corrupted shot inside of your request. Please try to remove the shot id for which the error first occurs from your `shotlist_request` and run DisruptionPy again. If the problem persists, please create an issue on the GitHub.


### Why does DisruptionPy continue to log `Processing result for shot: ***shot id***` after having retrieved all data? { .doc .doc-heading }
This is likeley a result of you `output_type_request` being slow. You should try to use a different `output_type_request` or batch the output of your current request.