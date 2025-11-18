"""manage port bundles for multiprocessing jobs

"""
import socket
from queue import Queue
import multiprocessing
import platform
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class PortBundlePool:
    """The Port Bundle Pool used by multiprocessing femag jobs
    """

    def __init__(self, host, port,
                 pool_size=3*multiprocessing.cpu_count(), port_bundle_size=3,
                 port_list=[]):
        self.host = host
        self.port = port
        self.pool_size = pool_size
        self.port_bundle_size = port_bundle_size
        self._pool = Queue(maxsize=pool_size)
        self.portList = port_list
        if not self.portList:
            self.initialize_pool()

    def initialize_pool(self):
        """Initialize the pool with pre-connected sockets."""
        logger.info(f'Initialize pool, size: {self.pool_size}, bundle_size: {self.port_bundle_size}')
        #if platform.system() == 'Windows':
        #    self.portList = [self.port + i * self.port_bundle_size
        #                     for i in range(self.pool_size)]
        #    logger.info(f'Initialize Done: {self.portList}')
        #    return
        # not windows
        while self._pool.qsize() < self.pool_size:
            logger.info(f'Initialize pool, size: {self._pool.qsize()}')
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            x = [self.test_socket(sock, self.port + i) for i in range(self.port_bundle_size)]
            if sum(x) == self.port_bundle_size:
                self._pool.put(self.port)
                self.portList.append(self.port)
                logger.debug(f'Found pool port bundle at: {self.port}')
            self.port += self.port_bundle_size
        logger.info(f'Initialize Done: {self.portList}')

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
