## Installation { .doc .doc-heading }

### How can I use DisruptionPy if I don't have access to the GitHub Repository? { .doc .doc-heading }
If you are unable to access the GitHub repository, you can manually install the package. Note that you may be installing an older version of DisruptionPy.

#### On CMod { .doc .doc-heading }
You can do this by running:
```bash
pip install /path/to/disruption-py
```

### Stuck on Poetry Install { .doc .doc-heading }
Poetry might hang instead of asking to unlock the keyring, which is a [known bug](https://github.com/python-poetry/poetry/issues/8623).
As a temporary workaround while the bug is addressed, please set:
```bash
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
```

## Usage { .doc .doc-heading }

### How do I debug issues with DisruptionPy? { .doc .doc-heading }
To better understand the issue that you are seeing in DisruptionPy you can enable more detailed logging by passing `log_settings` to the [`get_shots_data`][disruption_py.workflow.get_shots_data] call. For more information on customizing logging, please see [`LogSettings`][disruption_py.settings.LogSettings].

```python
import logging
from disruption_py.settings import LogSettings
from disruption_py.workflow import get_shots_data

get_shots_data(
    ...

    log_settings=LogSettings(
        log_to_console=True,
        log_file_path="path/to/log/file",
        file_log_level=logging.DEBUG,
        console_log_level=logging.DEBUG
    ),

    ...
)
```

### Why does DisruptionPy log `%MDSPLUS-E-ERROR` for all shots after a certain shot number? { .doc .doc-heading }
You may have a corrupted shot inside of your request. Please try to remove the shot id for which the error first occurs from your `shotlist_request` and run DisruptionPy again. If the problem persists, please create an issue on the GitHub.

### Why does DisruptionPy continue to log `Processing result for shot: ***shot id***` after having retrieved all data? { .doc .doc-heading }
This is likely a result of your [`output_setting`][disruption_py.settings.output_setting.OutputSetting] being slow. You should try to use a different `output_setting` or batch the output of your current request.