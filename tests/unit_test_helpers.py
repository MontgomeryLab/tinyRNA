"""
Helper functions for unit tests
"""
import subprocess
import threading
import hashlib
import signal
import psutil
import shlex
import sys
import io
import os


def get_dir_tree(root_path: str) -> dict:
    """Returns a nested dictionary representation of a given directory tree.

    Each directory is assigned its own sub-dictionary. Files within a directory
    are listed under the "files" set key within that directory's dictionary.
    This is very useful for checking correctness of operations which generate
    complex nested file outputs.

    Args:
        root_path: The relative or absolute path of the directory tree root.
            The base name (the lowest directory name defined in the path) will occupy
            the top level of the resulting dictionary.
    """

    root_dir = os.path.basename(os.path.realpath(root_path))
    dir_tree = {root_dir: {}}

    for top, dirs, files in os.walk(root_path):
        context = dir_tree[root_dir]
        if top != root_path:
            # This allows us to reference dir_tree dictionary keys relative to the current level and branch
            for level in top.replace(root_path + os.sep, '').split(os.sep):
                context = context[level]

        context['files'] = set(files)
        for folder in dirs: context[folder] = {}

    return dir_tree


def get_dir_checksum_tree(root_path: str) -> dict:
    """Returns a nested dictionary representation of a given directory tree with md5 file checksums

    Same as `get_dir_tree()` but each file in a directory's "files" set is represented as
    a tuple: (fileName, md5checksum). The md5 checksum is determined by the file's contents.
    This is useful for validating both structure AND content of complex nested file outputs.

    Args:
        root_path: The relative or absolute path of the directory tree root.
            The base name (the lowest directory name defined in the path) will occupy
            the top level of the resulting dictionary.
    """

    root_dir = os.path.basename(os.path.realpath(root_path))
    dir_tree = {root_dir: {}}

    for top, dirs, files in os.walk(root_path):
        context = dir_tree[root_dir]
        if top != root_path:
            # This allows us to reference dir_tree dictionary keys relative to the current level and branch
            for level in top.replace(root_path + os.sep, '').split(os.sep):
                context = context[level]

        context['files'] = set()
        for folder in dirs: context[folder] = {}
        for file in files:
            path = f"{top}{os.sep}{file}"
            with open(path, 'rb') as f:
                # Add (file, hash) tuple
                context['files'].add((file, hashlib.md5(f.read()).hexdigest()))

    return dir_tree


# Convenience function for resetting mocks between subtests
def reset_mocks(mock_open_f, mock_os, mock_stdout):
    mock_open_f.reset_mock()
    mock_os.reset_mock()
    mock_stdout.truncate(0)
    mock_stdout.seek(0)  # Avoids prepending a null string of previous buffer size


# The gzip module uses a buffered writer, so we need to reassemble the individual writes
def reassemble_gz_w(mock_calls):
    return b''.join([x[1][0] for x in mock_calls if x[0] == '().write'])


class ShellCapture:
    """Context manager for executing a shell command and reading the command's output

    Both blocking (execution of this script halts while the target completes) and
    non-blocking (execution of this script continues in parallel with the target) invocation
    is supported. Progress of the command can be measured using is_complete() and
    get_exit_status(). The entire child processes tree is properly terminated upon context
    exit to avoid creating zombie processes. STDOUT and STDERR are captured, meaning
    they are not printed to the console. Both STDOUT and STDERR can be retrieved for
    both blocking and non-blocking contexts.

    This object MUST be used as a context manager in a "with" block for proper behavior.
    However, its target command is executed with the () operator and you may interface
    with its attributes as you would any other object.

    Examples:
        with ShellCapture('for i in {1..5}; do echo "$i"; done') as shellEx:
            shellEx()
            print(f"The command output was {shellEx.get_stdout()}")

    Attributes:
        invocation: the command to be executed, with the appropriate blocking/non-blocking
            wrapper function to capture output streams
        blocking: a boolean variable to indicate which wrapper function was used
    """

    def __init__(self, command: str, blocking: bool = True) -> None:
        """Class constructor

        Args:
            command: The target shell command to be executed when this object is invoked with ()
            blocking: True for blocking execution, False for non-blocking execution
        """

        # Sanitize input to avoid any risk of shell injection, then tokenize
        tokenized_command = shlex.split(shlex.quote(command))
        # Assemble subprocess invocation
        self.blocking, self._result = blocking, None
        subprocess_method = subprocess.run if blocking else subprocess.Popen
        self.invocation = lambda: subprocess_method(
            tokenized_command,              # A list containing the command and each argument as separate elements
            stdout=subprocess.PIPE,         # Capture stdout
            stderr=subprocess.PIPE,         # Capture stderr
            universal_newlines=True,        # Assume utf-8 (non-binary) command output
            shell=True                      # Executes the command through a shell, allows for access to shell features
        )

    def __call__(self) -> None:
        """Allows this object to be called as though it were also a function"""

        self._result = self.invocation()

    def __enter__(self) -> 'ShellCapture':
        """Executed upon context entry.

        The default signal handler for SIGCHLD is saved so that it may be
        restored later. The default handler is replaced with a custom handler.
        """

        self._default_SIGCHLD_handler = signal.getsignal(signal.SIGCHLD)
        signal.signal(signal.SIGCHLD, self._handle_sigchld)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback) -> None:
        """Executed upon context exit or exception

        The entire child process tree is terminated from leaf to root, and
        the default SIGCHLD handler is restored.
        """

        for subproc in reversed(psutil.Process(os.getpid()).children(recursive=True)):
            try:
                subproc.terminate()
            except psutil.NoSuchProcess:
                pass  # If subprocess already exited, do nothing

        signal.signal(signal.SIGCHLD, self._default_SIGCHLD_handler)

    def get_stdout(self) -> str:
        """Retrieves the captured STDOUT of the target command"""

        if self._result is None: return ''
        elif self.blocking:
            return self._result.stdout  # subprocess.run -- stdout is a string object
        else:
            return self._result.communicate()[0]  # subprocess.Popen -- stdout is a stream object

    def get_stderr(self) -> str:
        """Retrieves the captured STDERR of the target command"""

        if self._result is None: return ''
        elif self.blocking:
            return self._result.stderr  # subprocess.run -- stderr is a string object
        else:
            return self._result.communicate()[1]  # subprocess.Popen -- stderr is a stream object

    def is_complete(self) -> bool:
        """Checks if the target command has finished execution"""

        if self.blocking:
            # A blocking call will have a non-None result only after complete execution
            return self._result is not None
        if not self.blocking and self._result is not None:
            # A non-blocking call may still be running. If so, .poll() will be None until complete
            return self._result.poll() is not None
        else:
            return False

    def get_exit_status(self) -> int:
        """Returns the target command's exit status if it has finished"""

        return self._result.returncode if self.is_complete() else None

    @staticmethod
    def _handle_sigchld(signum: int, stackframe: object) -> None:
        """Removes terminated child processes from the process table.

        This avoids creating zombie processes after termination (and the
        resource warnings that accompany)
        """

        try:
            os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            pass


class LambdaCapture:
    """Context manager for executing a Python function and reading the function's output

    Both blocking (execution of this script halts while the target completes) and
    non-blocking (execution of this script continues in parallel with the target) invocation
    is supported. Progress of the function can be measured using is_complete(). The entire
    child processes tree is properly terminated upon context exit to avoid creating zombie
    processes. Note: this behavior can be overridden by the target if the target uses the
    subprocess module or defines its own handler for SIGTERM or SIGCHLD. STDOUT and STDERR
    are captured, meaning they are not printed to the console. Note: this behavior can be
    overridden by the target. Both STDOUT and STDERR can be retrieved for both blocking and
    non-blocking contexts. The procedure for capturing stdout and stderr was informed by
    contextlib, which was not suitable for this application since daemon threads were exiting
    without proper context cleanup and this was leading to broken global stdout/stderr hooks.
    This implementation properly handles the transaction.

    This object MUST be used as a context manager in a "with" block for proper behavior.
    However, its target command is executed with the () operator and you may interface
    with its attributes as you would any other object.

    Examples:
        with LambdaCapture(lambda: print("oops", file=sys.stderr), blocking=False) as lambdaCap:
            lambdaCap()
            time.sleep(1)
            print(f"{lambdaCap.get_stderr()} was printed")

    Attributes:
        function_thread: The target function, wrapped in a thread which hooks stdout/stderr for capture
        blocking: A boolean variable to indicate whether the target function will be executed in
            serial or parallel
    """

    def __init__(self, function_call: callable, blocking: bool = True) -> None:
        """Class constructor

        Args:
            function_call: a lambda function to be executed when this object is invoked with ()
            blocking: True for blocking execution, False for non-blocking execution
        """

        self._stdout = io.StringIO()
        self._stderr = io.StringIO()
        self._old_targets = {}
        self.blocking = blocking

        # This is the thread target in which outputs are hooked before executing the provided function
        def redirect_target_output():
            # Save current stdout/stderr targets to restore on context exit
            self._old_targets = {stream: getattr(sys, stream) for stream in ["stdout", "stderr"]}
            # Set new targets to allow for capturing these streams
            setattr(sys, "stdout", self._stdout)
            setattr(sys, "stderr", self._stderr)
            # Finally, call target function
            function_call()

        # Ensure that a lambda function was provided to the constructor
        if isinstance(function_call, type(lambda:0)):
            # Wrap a thread around the capture function wrapped around the target function
            self.function_thread = threading.Thread(target=redirect_target_output, daemon=True)
        else:
            raise ValueError("Constructor requires a lambda function. Your argument was not.")

    def __call__(self) -> None:
        """Allows this object to be called as though it were also a function"""

        self.function_thread.start()
        if self.blocking: self.function_thread.join()

    def __enter__(self) -> 'LambdaCapture':
        """Executed upon context entry

        The default signal handler for SIGCHLD is saved so that it may be
        restored later. The default handler is replaced with a custom handler.
        """

        self._default_SIGCHLD_handler = signal.getsignal(signal.SIGCHLD)
        signal.signal(signal.SIGCHLD, self._handle_sigchld)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback) -> None:
        """Executed upon context exit or exception

        The entire child process tree is terminated from leaf to root, default
        STDOUT and STDERR behavior is restored, and the default SIGCHLD handler
        is restored.
        """

        for subproc in reversed(psutil.Process(os.getpid()).children(recursive=True)):
            try:
                subproc.terminate()
            except psutil.NoSuchProcess:
                pass  # If subprocess already exited, do nothing

        # Restore original stdout/stderr system attributes
        sys.stdout.flush()
        sys.stderr.flush()
        for stream, orig_target in self._old_targets.items():
            setattr(sys, stream, orig_target)

        # Restore SIGCHLD handler to default
        signal.signal(signal.SIGCHLD, self._default_SIGCHLD_handler)

    def get_stdout(self) -> str:
        """Retrieves the captured STDOUT of the target function"""

        if self.function_thread.ident is not None:
            return self._stdout.getvalue()
        else:
            return ''

    def get_stderr(self) -> str:
        """Retrieves the captured STDOUT of the target function"""

        if self.function_thread.ident is not None:
            return self._stderr.getvalue()
        else:
            return ''

    def is_complete(self) -> bool:
        """Returns true if the wrapper has been executed AND is no longer alive"""

        return self.function_thread.ident is not None and not self.function_thread.is_alive()

    @staticmethod
    def _handle_sigchld(signum: int, stackframe: object) -> None:
        """Removes terminated child processes from the process table.

        This avoids creating zombie processes after termination (and the
        resource warnings that accompany)
        """

        try:
            os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            pass
