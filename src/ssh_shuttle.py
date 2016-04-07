# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import paramiko
import os
import os.path as osp


class SSHShuttle(object):
    """
    Supports incremental or one-shot sending of files over to another
    server via SSH/SFTP.
    """

    def __init__(self, cfg):
        """
        Initialize via a configuration with keys:
          shuttle_sshkey, shuttle_remote_host, shuttle_remote_user and shuttle_remote_root
        
        :param cfg: dictionary with configuration keys
        """
        self.host = cfg['shuttle_remote_host']
        self.user = cfg['shuttle_remote_user']
        self.root = cfg['shuttle_remote_root']
        self.key = cfg['shuttle_ssh_key']

    
    def connect(self):
        """
        Initiate a connection with the remote server.
        """
        ssh = paramiko.SSHClient()
        ssh.load_system_host_keys()
        ssh.connect(self.host, username=self.user, key_filename=self.key)
        self.ssh = ssh
        self.sftp = ssh.open_sftp()

    
    def disconnect(self):
        """
        Close the connections gracefully.
        """
        self.sftp.close()
        self.ssh.close()

    
    def send_directory(self, local_dir, remote_path):
        """
        Send an entire directory of local files to the remote
        relative path.  Ensures remote directory exists.

        :param local_dir: the local directory to send
        :param remote_path: the remote path, appended to remote root if not absolute path
        """
        if remote_path[0] != '/':
           remote_path = osp.join(self.root, remote_path)

        self.chdir(remote_path, True)
        for file in os.listdir(local_dir):
            self.put(osp.join(local_dir, file), file)


    def put(self, local_path, remote_path):
        """
        Calls SFTP put to transfer the local file to the remote host.
        If the remote path is not absolute, the current directory on
        the remote host will be used.

        :param local_path: path to the local file
        :param remote_path: relative or absolute path to remote file
        """
        self.sftp.put(local_path, remote_path)


    def chdir(self, remote_dir, ensure_exists=False):
        """
        Change to the remote directory on target host.
        To mkdir the final part of the remote_dir if it does
        not exist, set ensure_exists.  Parent to remote_dir
        must always exist.

        :param remote_dir: the remote directory
        :param ensure_exists: make the final directory level if it does not exist
        """
        if ensure_exists:
            self.sftp.chdir(osp.dirname(remote_dir))
            try:
                self.sftp.mkdir(osp.basename(remote_dir))
            except IOError:
                # don't worry if directory already exists
                pass
        self.sftp.chdir(remote_dir)

