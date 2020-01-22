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
from builtins import str
from subprocess import check_call, check_output
import os
import os.path as osp
import re, sys
import json
import logging


class OutputCheckFailed(Exception):
    """
    Carries info about the failure of a subprocess execution.
    """
    pass


class Executor(object):
    """
    A class that handles execution of external processes for the system with
    make-like functionality.

    This class works for serially (not using MPI) and synchronously
    (not using a queuing system) launched executables.
    """

    def __init__(self, work_dir, exec_name):
        """
        Initialize the executor with the working directory, where the executable is located.

        :param work_dir: working directory, where the executable is located.
        :return:
        """
        self.work_dir = work_dir
        self.exec_name = exec_name

    def execute(self):
        """
        Execute the file in given working directory.

        The execute method first checks whether the expected output is present and
        if it indicates success, if yes, nothing is executed.  If not, the file is executed
        and its output is redirected: stdout goes to <exec_name>.stdout, stderr goes to <exec_name>.stderr
        Then the output is checked for markers indicating success, if found, function returns.
        Otherwise an exception is raised indicating the external process failed.

        :return: raises OutputCheckFailed if return code is non-zero
        """
        # first verify if we have already done our job
        try:
            self.check_output()
            return self
        except:
            pass

        logging.info('executing %s in %s' % (self.exec_name, self.work_dir))
        exec_name = self.exec_name
        stdout_file = open(osp.join(self.work_dir, exec_name + '.stdout'), 'w')
        stderr_file = open(osp.join(self.work_dir, exec_name + '.stderr'), 'w')
        logging.debug('EXECUTE: stdout %s' % stdout_file)
        logging.debug('EXECUTE: stderr %s' % stderr_file)
        check_call(exec_name, cwd=self.work_dir, stdout=stdout_file, stderr=stderr_file)
        return self

    def check_output(self):
        pass


class Geogrid(Executor):
    """
    Handles geogrid.exe execution.
    """

    def __init__(self, work_dir):
        """
        Only passes working directory.

        :param work_dir: working directory in which geogrid should be run
        :return:
        """
        super(Geogrid, self).__init__(work_dir, './geogrid.exe')

    def check_output(self):
        """
        Checks if the output file for geogrid contains the string "Successful completion of geogrid."

        :return: returns self on success, else raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')

        if not osp.exists(output_path):
            raise OutputCheckFailed("output file %s does not exist, cannot check output." % output_path)

        if 'Successful completion of geogrid.' not in open(output_path).read():
            print("Execution of %s was not successful." % self.exec_name)
            print("Examine %s for details." % output_path)
            raise OutputCheckFailed()

class Ungrib(Executor):
    """
    Handles ungrib.exe execution.
    """

    def __init__(self, work_dir):
        """
        Only passes working directory.

        :param work_dir:working directory in which ungrib should be run
        :return:
        """
        super(Ungrib, self).__init__(work_dir, './ungrib.exe')

    def check_output(self):
        """
        Checks if the output file for ungrib contains the string "Successful completion of ungrib."

        :return: raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of ungrib.' not in open(output_path).read():
            print(open(osp.join(self.work_dir, self.exec_name + '.stderr')).read())
            print('Execution of %s was not successful.' % self.exec_name)
            print('Examine %s for details.' % output_path)
            raise OutputCheckFailed()


class Metgrid(Executor):
    """
    Handles metgrid.exe execution.
    """

    def __init__(self, work_dir):
        """
        Only passes working directory.

        :param work_dir: working directory in which metgrid should be run
        :return:
        """
        super(Metgrid, self).__init__(work_dir, './metgrid.exe')

    def check_output(self):
        """
        Checks if the output file for metgrid contains the string "Successful completion of metgrid."

        :return: raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of metgrid.' not in open(output_path).read():
            print(open(osp.join(self.work_dir, self.exec_name + '.stderr')).read())
            print("Execution of %s was not successful." % self.exec_name)
            print("Examine %s for details." % output_path)
            raise OutputCheckFailed()


class Real(Executor):
    """
    Handles real.exe execution.
    """

    def __init__(self, work_dir):
        """
        Only passes working directory.

        :param work_dir: working directory in which real should be run
        :return:
        """
        super(Real, self).__init__(work_dir, './real.exe')

    def execute(self):
        """
        Execute real.exe in given working directory.

        This method is redefined here as real.exe needs special treatment.  DMPAR compilation of WRF
        causes real.exe to output stdout and stderr into rsl.out.0000 and rsl.error.0000 respectively.
        We don't redirect these (we can't) but we rename the output files after the fact.

        NOTE: on some machines it is OK to run real.exe from command line, but generally mpirun is required!

        :return: raises OutputCheckFailed if return code is non-zero
        """
        # first verify if we have already done our job
        try:
            self.check_output()
            return self
        except OutputCheckFailed:
            pass

        exec_name = self.exec_name
        wdir = self.work_dir
        check_call(exec_name, cwd=self.work_dir)
        os.rename(osp.join(wdir, "rsl.out.0000"), osp.join(wdir, "real.exe.stdout"))
        os.rename(osp.join(wdir, "rsl.error.0000"), osp.join(wdir, "real.exe.stderr"))

        return self

    def check_output(self):
        """
        Checks if the output file for real contains the string "SUCCESS COMPLETE REAL_EM INIT"

        :return: raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, "real.exe.stderr")

        if not osp.exists(output_path):
            raise OutputCheckFailed("output file %s does not exist, cannot check output." % output_path)

        if 'SUCCESS COMPLETE REAL_EM INIT' not in open(output_path).read():
            print("Execution of %s was not successful." % self.exec_name)
            print("Examine %s for details." % output_path)
            raise OutputCheckFailed()


class Submitter(object):
    """
    Class that abstract jobs that must be submitted to a queue manager.
    """

    def __init__(self, work_dir, qsys_id):
        """
        Initialize with working directory and queue manager id.

        :param work_dir: where the job is to be run
        :param qsys_id: the id of the queue system/template to use (points to etc/clusters.json)
        """
        self.qsys_infos = json.load(open('etc/clusters.json'))
        self.work_dir = work_dir
        if qsys_id not in self.qsys_infos:
            raise ValueError('Invalid queue system, must be one of %s' % repr(list(self.qsys_infos.keys())))
        self.qsys_id = qsys_id


    def submit(self, task_id, exec_path, nodes, ppn, wall_time_hrs):
        """
        Build a job script and submit it to the queue manager.

        :param task_id: the name of the task in the queue
        :param exec_path: the path to the parallel executable
        :param nodes: number of nodes to request for parallel job
        :param ppn: processors per nodes to request
        :param wall_time_hrs: wall time to request for job
        :return: job number string to be used in further queue manager commands
        """
        qsys = self.qsys_infos[self.qsys_id]
        qsub = qsys['qsub_cmd']
        script_tmpl = open(qsys['qsub_script']).read()

        args = {'nodes': nodes, 'ppn': ppn, 'wall_time_hrs': wall_time_hrs,
                'exec_path': exec_path, 'cwd': self.work_dir, 'task_id': task_id, 'np': nodes * ppn}

        script_path = osp.join(self.work_dir, task_id + ".sub")
        with open(script_path, 'w') as f:
            f.write(script_tmpl % args)

        logging.info('submitting to batch queue system %s' % self.qsys_id)
        logging.info('%s %s' % (qsub, script_path))

        ret = check_output([qsub, script_path], cwd=self.work_dir).decode(sys.stdout.encoding)
        logging.info(ret)
        job_num = ret.split(' ')[qsys['qsub_job_num_index']].rstrip()
        logging.info('job number %s submitted' % job_num)
        return job_num


class WRF(Submitter):
    """
    Handles job submission for wrf.exe into a job manager.
    """

    def __init__(self, work_dir, qman):
        """
        Initialize with working directory and queue manager to interface with.

        :param work_dir: working directory
        :param qman: the queue manager code (sge, ...)
        """
        super(WRF, self).__init__(work_dir, qman)

    def submit(self, task_id, nodes, ppn, wall_time_hrs):
        """
        Submits a WRF job with the given parameters via queue manager setup during construction.

        Fills out executable name and lets submitter do the heavy lifting.

        :param task_id: the name of the task in the queue
        :param nodes: number of nodes to request for parallel job
        :param ppn: processors per nodes to request
        :param wall_time_hrs: wall time to request for job
        :return: job number string to be used in further queue manager commands
        """
        ret = super(WRF, self).submit(task_id, "./wrf.exe", nodes, ppn, wall_time_hrs)
        return ret

