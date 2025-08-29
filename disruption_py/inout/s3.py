import threading
from loguru import logger
import s3fs

from disruption_py.machine.tokamak import Tokamak
from disruption_py.config import config


class S3Connection:
    """
    Class for connecting to S3 bucket.
    """

    def __init__(self, endpoint_url: str):
        self.conn = None
        if endpoint_url is None:
            return
        logger.debug(
            "PID #{pid} | Connecting to S3 bucket: {server}",
            server=endpoint_url,
            pid=threading.get_native_id(),
        )
        # pylint: disable=no-member
        self.conn = s3fs.S3FileSystem(anon=True, endpoint_url=endpoint_url)

    @classmethod
    def from_config(cls, tokamak: Tokamak):
        """
        Create instance of the MDS connection based on the connection string
        from the configuration.
        """
        return S3Connection(config(tokamak).inout.s3.endpoint_url)

    def get_shot_connection(self, shot_id: int):
        """Get S3 Connection wrapper for individual shot."""
        return self

    def cleanup(self):
        pass
