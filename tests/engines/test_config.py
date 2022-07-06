import femagtools.config
import os.path
import platform


def test_femag():
    if platform.system() == "Windows":
        expected = "wfemagw64"
    else:
        expected = "xfemag64"

    try:
        msg = femagtools.config.get_executable()
    except Exception as e:
        msg = str(e)
    assert msg.find(expected) > -1


def test_engine_config():
    default_config = {
        'ENGINE': 'amazon',
        'INSTANCE_TYPE': 't2.micro'}
    config = femagtools.config.Config(default_config)
    ini_file = os.path.join(os.path.dirname(__file__), 'config.ini')
    config.from_ini_file(ini_file)
    assert config['KEY_NAME'] == 'test-key-pair-eucentral-1'
