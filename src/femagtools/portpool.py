"""manage port bundles for multiprocessing jobs

"""
import socket
from queue import Queue
import multiprocessing
import platform
import logging

logger = logging.getLogger(__name__)

class PortBundlePool:
    """The Port Bundle Pool used by multiprocessing femag jobs
    """

    def __init__(self, host, port,
                 pool_size=3*multiprocessing.cpu_count(), port_bundle_size=3,
                 port_list=[],
                 notify=None):
        self.host = host
        self.port = port
        self.pool_size = pool_size
        self.port_bundle_size = port_bundle_size
        self._pool = Queue(maxsize=pool_size) # not really used
        self.portList = port_list
        self.notify = notify
        logger.info(f'PortBundlePool PortList: {self.portList}')
        if not self.portList:
            self.initialize_pool()

    def initialize_pool(self):
        """Initialize the pool with pre-connected sockets."""
        logger.info(f'Initialize pool, size: {self.pool_size}, bundle_size: {self.port_bundle_size}')
        # not windows
        logger.info(f'Initialize Multiprocessing Port Pool')
        if self.notify:
            publish_main = f'<progress>Initialize Multiprocessing Port Pool'
            self.notify(['femag_log', publish_main])
        while self._pool.qsize() < self.pool_size:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            x = [self.test_socket(sock, self.port + i) for i in range(self.port_bundle_size)]
            if sum(x) == self.port_bundle_size:
                self._pool.put(self.port)
                self.portList.append(self.port)
                logger.debug(f'Found pool port bundle at: {self.port}')
                if self.notify:
                    self.notify(['femag_log', f'{publish_main}<br/>Found port bundle at: {self.port}'
                                 f'<br/>Port {self._pool.qsize()} of {self.pool_size}'])
            self.port += self.port_bundle_size
        logger.info(f'Initialize Done: {self.portList}')
        if self.notify:
            self.notify(['femag_log', f'Initialize Multiprocessing Port Pool: {self.portList}'])

    def get_port_list(self):
        """Returns a list of free/not used socket ports"""
        return self.portList

    def test_socket(self, sock, port):
        """Test a socket"""
        try:
            sock.connect((self.host, port))
            logger.debug(f'Port: {port} is open')
            return False
        except:
            pass
        return True

# Example usage
if __name__ == "__main__":
    # port pool
    pool = PortBundlePool("localhost", 8880, pool_size=3*multiprocessing.cpu_count(), port_bundle_size=3)
    print(f'pool port list: {pool.get_port_list()}')
