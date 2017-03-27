from ssh_shuttle import SSHShuttle
from forecast import load_sys_cfg
import sys

cfg = load_sys_cfg()
s = SSHShuttle(cfg)
s.connect()
command = ' '.join(sys.argv[1:])
s.simple_command(command)
