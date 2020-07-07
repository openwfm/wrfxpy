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


from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
import paramiko
import os
import os.path as osp
import logging
import json
import glob
import sys 
import pprint
import fcntl
import errno
import collections
from utils import load_sys_cfg



class SSHShuttle(object):
    """
    Supports incremental or one-shot sending of files over to another
    server via SSH/SFTP.
    """

    def __init__(self, cfg):
        """
        Initialize via a configuration with keys:
          shuttle_sshkey, shuttle_remote_host, shuttle_remote_user and shuttle_remote_root
          plus anything else that shuttle will need from cfg
        
        :param cfg: dictionary with configuration keys
        """
        self.host = cfg['shuttle_remote_host']
        self.user = cfg['shuttle_remote_user']
        self.root = cfg['shuttle_remote_root']
        self.key = cfg['shuttle_ssh_key']
        self.lock_path = cfg['shuttle_lock_path']
        self.lock_file = None
        self.workspace_path = cfg['workspace_path']
        self.cat_local_path = osp.join(self.workspace_path, 'catalog.json')
        self.connected=False

    
    def connect(self):
        """
        Initiate a connection with the remote server.
        """
        logging.info('SHUTTLE connecting to remote host %s' % self.host)
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.load_system_host_keys()
        #ssh.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))
        ssh.connect(self.host, username=self.user, key_filename=self.key)
        self.ssh = ssh
        self.sftp = ssh.open_sftp()
        # change to the remote root immediately
        logging.info('SHUTTLE changing to remote root %s' % self.root)
        self.sftp.chdir(self.root)
        self.connected = True

    
    def disconnect(self):
        """
        Close the connections gracefully.
        """
        logging.info('SHUTTLE operations completed, closing connection.')
        self.sftp.close()
        self.ssh.close()
        self.connected = False

    
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
        if remote_path[0] != '/':
           remote_path = osp.join(self.root, remote_path)
        self.sftp.put(local_path, remote_path)

    
    def get(self, remote_path, local_path):
        """
        Uses an SFTP get operation to transfer the remote file to
        the local filesystem.  The remote path can be absolute or 
        relative to the stored root.

        :param remote_path: remote path
        :param local_path: local path
        """
        if remote_path[0] != '/':
           remote_path = osp.join(self.root, remote_path)
        try:
            self.sftp.get(remote_path, local_path)
        except IOError as e:
            logging.error('ssh_shuttle: sftp.get failed: %s' % e.strerror)
            logging.error('remote_path: %s' % remote_path)
            logging.error('local_path: %s' % local_path)
            exit(1)


    def chdir(self, remote_dir, ensure_exists=False):
        """
        Change to the remote directory on target host.
        To mkdir the final part of the remote_dir if it does
        not exist, set ensure_exists.  Parent to remote_dir
        must always exist.

        :param remote_dir: the remote directory
        :param ensure_exists: make the final directory level if it does not exist
        """
        if remote_dir[0] != '/':
           remote_dir = osp.join(self.root, remote_dir)
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
        logging.info('SHUTTLE removing directory ' + abs_path)
        self.simple_command('/bin/rm -rf ' + abs_path)
        #files = self.sftp.listdir(abs_path) # was taking too long
        #for f in files:
        #    self.sftp.remove(osp.join(abs_path, f))
        #self.sftp.rmdir(abs_path)

    def simple_command(self, command):
        """
        Execute command on the remote host with no frills.

        :param command: the command string to be executed
        """
        stdin, stdout, stderr = self.ssh.exec_command(command)
        stdin.flush()
        print(stdout.read().decode())
        print(stderr.read().decode())

    def retrieve_catalog(self):
        """ 
        Retrieve the catalog from the remote visualization server
        and update the local copy
        NOTE: connect and acquire lock before catalog operations
        
        """
        logging.info('SHUTTLE retrieving catalog file.')
        self.simple_command('wrfxweb/join_catalog.sh')
        self.get('catalog.json', self.cat_local_path)
        cat = json.load(open(self.cat_local_path))
        logging.info('SHUTTLE retrieve complete.')
        return cat 

def send_product_to_server(cfg, local_dir, remote_dir, sim_name, manifest_filename, description = None, exclude_files = []):
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
    
    logging.info('SHUTTLE send_product_to_server')
    # logging.info('SHUTTLE local directory    %s' % local_dir)
    logging.info('SHUTTLE remote directory   %s' % remote_dir)
    logging.info('SHUTTLE simulation name    %s' % sim_name)
    logging.info('SHUTTLE manifest file name %s' % manifest_filename)
    logging.info('SHUTTLE description        %s' % description)
    logging.debug('SHUTTLE configuration:\n%s' % pprint.pformat(cfg,indent=4))
   
    s = SSHShuttle(cfg)
    s.connect()

    logging.info('SHUTTLE root directory     %s' % s.root)
    logging.info('SHUTTLE sending local directory %s to remote host' % local_dir)
    sent_files = s.send_directory(local_dir, remote_dir, exclude_files)

    # identify the start/end UTC time (all domains may not have the same simulation extent)
    # if more than one manifest file match, take the first one
    manifest_file = glob.glob(osp.join(local_dir,manifest_filename)) 
    if len(manifest_file) == 0:
        logging.warning('SHUTTLE did not find manifest file %s, catalog update postponed' % manifest_filename)
    else:
        manifest_file = manifest_file[0]
        logging.info('SHUTTLE loading manifest filename given %s loading %s' % (manifest_filename, manifest_file))
        mf = json.load(open(manifest_file))
        #print 'manifest:', json.dumps(mf, indent=4, separators=(',', ': '))
        logging.debug('manifest %s' % mf)
        k=[]
        for i in mf.keys():
            dom = mf[i]
            k = k + list(dom.keys())
        times = sorted(k)
        logging.info('SHUTTLE detected local start/end UTC times as %s - %s' % (times[0], times[-1]))
    
        # retrieve the catalog & update it
        logging.info('SHUTTLE updating local catalog file on remote host')
        local_cat = { sim_name : { 'manifest_path' : '%s/%s' % (remote_dir, osp.basename(manifest_file)),
                          'description' : description if description is not None else sim_name,
                          'from_utc' : times[0],
                          'to_utc' : times[-1] }
                    }
        local_cat_path = osp.join(local_dir,'catalog.json')
        json.dump(local_cat, open(local_cat_path, 'w'), indent=1, separators=(',',':'))
        remote_cat_path = remote_dir + '/catalog.json'
        s.put(local_cat_path, remote_cat_path)
        # s.simple_command('ls -l %s' % osp.join(s.root,remote_cat_path))
        s.simple_command('wrfxweb/join_catalog.sh')

    s.disconnect()

    logging.info('SHUTTLE sent %d files to visualization server.'  % len(sent_files))
    return sent_files

def ssh_command(command):
    cfg = load_sys_cfg()
    s = SSHShuttle(cfg)
    s.connect()
    s.simple_command(command)
    s.disconnect()

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print(('usage: %s <local-dir> <remote-relative-dir> <sim-name>' % sys.argv[0]))
        sys.exit(1)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    local_dir = sys.argv[1]
    remote_dir = sys.argv[2]
    sim_name = sys.argv[3]

    cfg = load_sys_cfg()

    send_product_to_server(cfg, local_dir, remote_dir, sim_name, None, [])

