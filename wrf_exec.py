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


class OutputCheckFailed(Exception):
    """
    Carries info about the failure of a subprocess execution.
    """
    pass


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

        :return: raises OutputCheckFailed if return code is non-zero
        """
        # first verify if we have already done our job
        try:
            self.check_output()
            return self
        except:
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

        :return: returns self on success, else raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')

        if not osp.exists(output_path):
            raise OutputCheckFailed("output file %s does not exist, cannot check output." % output_path)

        if 'Successful completion of geogrid.' not in open(output_path).read():
            raise OutputCheckFailed(
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

        :return: raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of ungrib.' not in open(output_path).read():
            raise OutputCheckFailed(
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

        :return: raises OutputCheckFailed
        """
        output_path = osp.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of metgrid.' not in open(output_path).read():
            raise OutputCheckFailed(
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
            raise OutputCheckFailed(
                    "Execution of real.exe was not successful, examine %s for details." % output_path)


class Submitter(object):
    """
    Class that abstract jobs that must be submitted to a queue manager.
    """

    def __init__(self, work_dir, qman):
        """
        Initialize with working directory and queue manager id.

        :param work_dir: where the job is to be run
        :param qman: id of queue manager, affects job submission script construction
        """
        self.work_dir = work_dir
        if qman not in self.qman_infos:
            raise ValueError('Invalid queue manager, must be one of %s' % repr(self.queue_managers.keys()))
        self.qman = qman

    def submit(self, task_id, exec_path, nodes, ppn, wall_time_hrs):
        """
        Build a job script and submit it to the queue manager.

        :param task_id: the name of the task in the queue
        :param exec_path: the path to the parallel executable
        :param nodes: number of nodes to request for parallel job
        :param ppn: processors per nodes to request
        :param wall_time_hrs: wall time to request for job
        """
        qman_info = self.qman_infos[self.qman]
        qsub = qman_info['qsub_cmd']
        script_tmpl = qman_info['qsub_script']
        job_code_f = qman_info['job_code_f']

        args = {'nodes': nodes, 'ppn': ppn, 'wall_time_hrs': wall_time_hrs,
                'exec_path': exec_path, 'cwd': self.work_dir, 'task_id': task_id, 'np': nodes * ppn}

        script_path = osp.join(self.work_dir, task_id + ".sub")
        with open(script_path, 'w') as f:
            f.write(script_tmpl % args)

        ret = check_call([qsub, script_path], cwd=self.work_dir)
        try:
            self.job_code = job_code_f(ret)
        except ValueError as e:
            raise CalledProcessException('Failed to capture job code from submit with error %s' % e)

    qman_infos = {
        'sge': {
            'qsub_cmd': 'qsub',
            'job_code_f': lambda x: int(x.strip()),
            'qsub_script': '#$ -S /bin/bash\n' '#$ -N %(task_id)s\n' '#$ -wd %(cwd)s\n'
                           '#$ -l h_rt=%(wall_time_hrs)02d:00:00\n' '#$ -pe mpich %(np)d\n'
                           'mpirun_rsh -np %(np)d -hostfile $TMPDIR/machines %(exec_path)s\n'
        }
    }


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
        """
        super(WRF, self).submit(task_id, "./wrf.exe", nodes, ppn, wall_time_hrs)
