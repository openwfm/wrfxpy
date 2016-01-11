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

from subprocess import check_call 
import os
import os.path as osp

class ExecFailedError(Exception):
    """
    Carries info about the failure of a subprocess execution.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Executor(object):
    """
    A class that handles execution of external processes for the system.

    This class works for executables where we first run them and then check the output
    for indications of success.
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

        Output redirection: stdout goes to <exec_name>.stdout, stderr goes to <exec_name>.stderr

        :return: raises ExecFailedError if return code is non-zero
        """
        # first verify if we have already done our job
        try:
          self.check_output()
          return self
        except ExecFailedError:
          pass

        exec_name = self.exec_name
        stdout_file = open(osp.join(self.work_dir, exec_name + '.stdout'), 'w')
        stderr_file = open(osp.join(self.work_dir, exec_name + '.stderr'), 'w')
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

        :return: returns self on success, else raises ExecFailedError
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of geogrid.' not in open(output_path).read():
            raise ExecFailedError(
                "Execution of geogrid.exe was not successfull, examine %s for details." % output_path)


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

        :return: raises ExecFailedError
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of ungrib.' not in open(output_path).read():
            raise ExecFailedError(
                "Execution of ungrib.exe was not successful, examine %s for details." % output_path)


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

        :return: raises ExecFailedError
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of metgrid.' not in open(output_path).read():
            raise ExecFailedError(
                "Execution of ungrib.exe was not successful, examine %s for details." % output_path)


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

        :return: raises ExecFailedError if return code is non-zero
        """
        # first verify if we have already done our job
        try:
          self.check_output()
          return self
        except ExecFailedError:
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

        :return: raises ExecFailedError
        """
        output_path = osp.join(self.work_dir, "real.exe.stderr")
        if 'SUCCESS COMPLETE REAL_EM INIT' not in open(output_path).read():
            raise ExecFailedError(
                "Execution of real.exe was not successful, examine %s for details." % output_path)



class WRF(Executor):
    """
    Handles job submission for wrf.exe into a job manager.
    """

    def __init__(self, work_dir, qman):
        """
        Initialize with working directory and queue manager to interface with.

        :param work_dir: working directory
        :param qman: the queue manager code (sge, ...)
        """
        super(WRF, self).__init__(work_dir)
        if qman not in self.queue_managers:
            raise ValueError('Invalid queue manager, must be one of %s' % repr(self.queue_managers.keys()))
        self.qman = qman


    def submit(self, task_id, nodes, ppn, wall_time_hrs):
        """
        Submit the job into the queue manager with given parameters.
        """
        args = { 'nodes': nodes, 'ppn': ppn, 'wall_time_hrs': wall_time_hrs,
                 'cwd': self.work_dir, 'task_id': task_id, 'np': nodes * ppn }
                 


    qman_scripts = { 'sge' : """
#$ -S /bin/bash
#$ -N %(task_id)s
#$ -wd %(cwd)s
#$ -l h_rt=%(wall_time_hrs)02d:00:00
#$ -pe mpich %(np)d
mpirun_rsh -np %(np)d -hostfile $TMPDIR/machines ./wrf.exe
"""
 }

