from __future__ import absolute_import
from ssh_shuttle import ssh_command 

command = ' '.join(sys.argv[1:])
ssh_command(command)
