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
import logging
import json
import glob


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
        # change to the remote root immediately
        self.sftp.chdir(self.root)

    
    def disconnect(self):
        """
        Close the connections gracefully.
        """
        self.sftp.close()
        self.ssh.close()

    
    def send_directory(self, local_dir, remote_path, exclude_files = []):
        """
        Send an entire directory of local files to the remote
        relative path.  Ensures remote directory exists.

        :param local_dir: the local directory to send
        :param remote_path: the remote path, appended to remote root if not absolute path
        :param exclude_files: files to be excluded from the directory, default is []
        :return: the filenames that were put via SFTP
        """
        if remote_path[0] != '/':
           remote_path = osp.join(self.root, remote_path)

        prev_cwd = self.sftp.getcwd()
        self.chdir(remote_path, True)

        exclude_set = set(exclude_files)
        sent_files = []
        for filename in os.listdir(local_dir):
            if filename not in exclude_set:
	            self.sftp.put(osp.join(local_dir, filename), filename)
	            sent_files.append(filename)

        # restore previous working directory
        self.chdir(prev_cwd)
        return sent_files


    def put(self, local_path, remote_path):
        """
        Calls SFTP put to transfer the local file to the remote host.

        :param local_path: path to the local file
        :param remote_path: relative or absolute path to remote file
        """
        self.sftp.put(local_path, remote_path)

    
    def get(self, remote_path, local_path):
        """
        Uses an SFTP get operation to transfer the remote file to
        the local filesystem.  The remote path can be absolute or 
        relative to the stored root.

        :param remote_path: remote path
        :param local_path: local path
        """
        self.sftp.get(remote_path, local_path)


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


    def rmdir(self, remote_dir):
        """
        Remove a remote directory.

        :param remote_dir: the directory to remove
        """
        abs_path = osp.join(self.root, remote_dir)
        files = self.sftp.listdir(abs_path)
        for f in files:
            self.sftp.remove(osp.join(abs_path, f))
        self.sftp.rmdir(abs_path)


def send_product_to_server(cfg, local_dir, remote_dir, sim_name, description = None, exclude_files = []):
    """
    Executes all steps required to send a local product directory to a remote visualization
    server for display.

    :para cfg: the configuration file for the shuttle
    :param local_dir: the local product directory
    :param remote_dir: the remote directory (relative to remote root set in conf.json)
    :param sim_name: the simulation name and id
    :param description: the description to be placed into the catalog file
    :param exclude_files: filenames that are not skipped during the upload, default is []
    :return: a list of the files that was sent
    """
    s = SSHShuttle(cfg)
    logging.info('SHUTTLE connecting to remote host %s' % s.host)
    s.connect()

    # identify the catalog file
    manifest_pattern = osp.join(local_dir, '*.json')
    logging.info('looking for local manifest file %s' % manifest_pattern)
    manifest_file = glob.glob(osp.join(local_dir, '*.json'))[0]
    logging.info('SHUTTLE found local manifest file %s' % manifest_file)

    logging.info('SHUTTLE sending local direcory %s to remote host' % local_dir)
    sent_files = s.send_directory(local_dir, remote_dir, exclude_files)

    # identify the start/end UTC time (all domains have same simulation extent)
    mf = json.load(open(manifest_file))
    dom = mf[mf.keys()[0]]
    times = sorted(dom.keys())
    logging.info('SHUTTLE detected local start/end UTC times as %s - %s' % (times[0], times[-1]))

    # retrieve the catalog & update it
    logging.info('SHUTTLE updating catalog file on remote host')
    cat_local = osp.join(cfg['workspace_path'], 'catalog.json')
    s.get('catalog.json', cat_local)

    cat = json.load(open(cat_local))
    cat[sim_name] = { 'manifest_path' : '%s/%s' % (remote_dir, osp.basename(manifest_file)),
                      'description' : description if description is not None else sim_name,
                      'from_utc' : times[0],
                      'to_utc' : times[-1] }
    json.dump(cat, open(cat_local, 'w'), indent=4, separators=(',', ': '))
    s.put(cat_local, 'catalog.json')

    logging.info('SHUTTLE operations completed, closing connection.')

    s.disconnect()
    return sent_files


if __name__ == '__main__':
    import sys 

    if len(sys.argv) != 4:
        print('usage: %s <local-dir> <remote-relative-dir> <sim-name>' % sys.argv[0])
        sys.exit(1)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    local_dir = sys.argv[1]
    remote_dir = sys.argv[2]
    sim_name = sys.argv[3]

    cfg = json.load(open('etc/conf.json'))

    send_product_to_server(cfg, local_dir, remote_dir, sim_name, None, [])
 
