import femagtools.config
import os.path


def test_config():
    default_config = {
        'ENGINE': 'amazon',
        'INSTANCE_TYPE': 't2.micro'}
    config = femagtools.config.Config(default_config)
    ini_file = os.path.join(os.path.dirname(__file__), 'config.ini')
    config.from_ini_file(ini_file)
    assert config['KEY_NAME'] == 'test-key-pair-eucentral-1'
