#!/usr/bin/env python
#
import unittest
import femagtools
try:
    import mock
except ImportError:
    import unittest.mock as mock
import os


class PropulationTest(unittest.TestCase):

    @mock.patch('femagtools.amazon.Engine._create_amazon_resource')
    def setUp(self, create_amazone_resource):
        create_amazone_resource = mock.Mock()
        self.engine = femagtools.amazon.Engine()

    def tearDown(self):
        self.engine = None

    def test_read_cloud_ini(self):
        testPath = os.path.join(os.path.split(__file__)[0], 'data')
        if len(testPath) == 0:
            testPath = os.path.join(os.path.abspath('.'), 'data')

        result = """Test Cloud init
export ACL=authenticated-read
export CLOUD_INIT={}/cloud_init.txt
export COMPANY_NAME=femag
export ENGINE=amazon
export FINISH_TASK_FILENAME=exit_code
export IMAGE_ID=ami-b0cc23df
export INSTANCE_TYPE=t2.micro
export SERVER_LOCATION=eu-central-1
export BUCKET_NAME=1
""".format(testPath)

        self.engine.config['CLOUD_INIT'] = "{}/cloud_init.txt".format(testPath)
        user_data = self.engine._read_cloud_init('1')
        self.assertEqual(user_data, result)

    # This test does not work, cause mock can not mock attributes in a class which are
    # not defined in the init method
    @mock.patch('femagtools.job.CloudJob')
    @mock.patch('femagtools.job.Task')
    def test_get_status_code(self, MockCloudJob, MockCloudTask):
        filename = "../data/exit_code"
        from femagtools.job import CloudTask

        # MockCloudJob.add_task = Mock()
        # MockCloudTask.directory.__get__ = Mock(return_value='/tmp')
        # MockCloudJob.tasks.return_value = [MockCloudTask, MockCloudTask, MockCloudTask]

        m = mock.MagicMock()
        mock.directory = mock.Mock(return_value='/tmp')
        self.engine.create_job('')

        self.engine.job.tasks.append(m(2, "/tmp/"))

        # This is the test
        # self.engine._get_status_code('exit_code')


if __name__ == '__main__':
    unittest.main()
