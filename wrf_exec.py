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

from subprocess import check_call, CalledProcessError
import os.path as pth


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

        :return: raises CalledProcessError if return code is non-zero
        """
        exec_name = self.exec_name
        stdout_file = open(pth.join(self.work_dir, exec_name + '.stdout'), 'w')
        stderr_file = open(pth.join(self.work_dir, exec_name + '.stderr'), 'w')
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

        :return: returns self on success, else raises CalledProcessError
        """
        output_path = pth.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of geogrid.' not in open(output_path).read():
            raise CalledProcessError(
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

        :return: raises CalledProcessError
        """
        output_path = pth.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of ungrib.' not in open(output_path).read():
            raise CalledProcessError(
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
        super(Ungrib, self).__init__(work_dir, './metgrid.exe')


    def check_output(self):
        """
        Checks if the output file for metgrid contains the string "Successful completion of metgrid."

        :return: raises CalledProcessError
        """
        output_path = pth.join(self.work_dir, self.exec_name + '.stdout')
        if 'Successful completion of metgrid.' not in open(output_path).read():
            raise CalledProcessError(
                "Execution of ungrib.exe was not successful, examine %s for details." % output_path)


