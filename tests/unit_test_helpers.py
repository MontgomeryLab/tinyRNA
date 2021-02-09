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

"""
Returns a nested dictionary representation of a given directory tree.
Each directory is assigned its own sub-dictionary. Files within a
directory are listed under the "files" array key within that directory's
dictionary. This is very useful for checking correctness of operations
which generate complex nested file outputs.
"""


def get_dir_tree(root_path):
    root_dir = os.path.basename(root_path)
    dir_tree = {root_dir: {}}

    for top, dirs, files in os.walk(root_path):
        context = dir_tree[root_dir]
        if top != root_path:
            # This allows us to reference dir_tree dictionary keys relative to the current level and branch
            for level in top.replace(root_path + os.sep, '').split(os.sep):
                context = context[level]

        context['files'] = files
        for folder in dirs: context[folder] = {}

    return dir_tree


"""
Returns a nested dictionary representation of a given directory tree, as above.
However, this function also tags each file with an md5 checksum of its contents.
This is useful for validating both structure AND content of complex nested
file outputs.
"""


def get_dir_checksum_tree(root_path):
    root_dir = os.path.basename(root_path)
    dir_tree = {root_dir: {}}

    for top, dirs, files in os.walk(root_path):
        context = dir_tree[root_dir]
        if top != root_path:
            # This allows us to reference dir_tree dictionary keys relative to the current level and branch
            for level in top.replace(root_path + os.sep, '').split(os.sep):
                context = context[level]

        context['files'] = []
        for folder in dirs: context[folder] = {}
        for file in files:
            path = f"{top}/{file}"
            with open(path, 'rb') as f:
                # Add (file, hash) tuple
                context['files'].append((file, hashlib.md5(f.read()).hexdigest()))

    return dir_tree


"""
Context manager for executing the given post-install aquatx command and
reading the command's stdout/stderr, then properly terminating any remaining
children to avoid zombie processes. It's not as dark as it sounds.
"""


class ShellCapture:
    # Constructor
    def __init__(self, command, blocking=True):
        # Sanitize input to avoid any risk of shell injection, then tokenize
        tokenized_command = shlex.split(shlex.quote(command))
        # Assemble subprocess invocation
        self.blocking, self.result = blocking, None
        subprocess_method = subprocess.run if blocking else subprocess.Popen
        self.invocation = lambda: subprocess_method(
            tokenized_command,              # A list containing the command and each argument as separate elements
            stdout=subprocess.PIPE,         # Capture stdout
            stderr=subprocess.PIPE,         # Capture stderr
            universal_newlines=True,        # Assume utf-8 (non-binary) command output
            shell=True                      # Executes the command through a shell, allows for access to shell features
        )

    # Allows this object to be called as though it were also a function
    def __call__(self):
        # Set handler for child process termination (avoids zombie children)
        signal.signal(signal.SIGCHLD, self.__handle_sigchld)
        self.result = self.invocation()

    # Executed upon context entry
    def __enter__(self):
        # Save the default SIGCHLD handler to restore later
        self.default_SIGCHLD_handler = signal.getsignal(signal.SIGCHLD)
        return self

    # Executed upon context exit or exception
    def __exit__(self, exception_type, exception_value, exception_traceback):
        # Remove child processes starting with lowest descendents first
        for subproc in reversed(psutil.Process(os.getpid()).children(recursive=True)):
            subproc.terminate()

        # Restore SIGCHLD handler to default
        signal.signal(signal.SIGCHLD, self.default_SIGCHLD_handler)

    def get_stdout(self):
        if self.result is None: return ''
        elif self.blocking:  # subprocess.run -- stdout is a string object
            return self.result.stdout
        else:  # subprocess.Popen -- stdout is a stream object
            return self.result.communicate()[0]

    def get_stderr(self):
        if self.result is None: return ''
        elif self.blocking:  # subprocess.run -- stderr is a string object
            return self.result.stderr
        else:  # subprocess.Popen -- stderr is a stream object
            return self.result.communicate()[1]

    def is_complete(self):
        if self.blocking:
            # A blocking call will have a non-None result only after complete execution
            return self.result is not None
        if not self.blocking and self.result is not None:
            # A non-blocking call may still be running. If so, .poll() will be None until complete
            return self.result.poll() is not None
        else:
            return False

    def get_exit_status(self):
        return self.result.returncode if self.is_complete() else None

    @staticmethod
    def __handle_sigchld(signum, stackframe):
        # Removes terminated child processes from the process table.
        # This avoids creating zombie processes after termination.
        try:
            os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            pass


"""
Context manager for executing the given pre-install Python function
in a new thread and reading the resulting stdout/stderr. The procedure
for capturing stdout and stderr was informed by contextlib, which was not
suitable for this application since daemon threads were exiting without
proper context cleanup and this was leading to broken global stdout/stderr
hooks. This implementation properly handles the transaction.
"""


class LambdaCapture:

    # Constructor takes a lambda function to be executed
    def __init__(self, function_call, blocking=True):
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

    # Allows this object to be called as though it were also a function
    def __call__(self):
        # Set handler for child process termination (avoids zombie children)
        signal.signal(signal.SIGCHLD, self.__handle_sigchld)
        # Execute the target function
        self.function_thread.start()
        # Wait for the thread to finish before returning from this call
        if self.blocking: self.function_thread.join()

    # Executed upon context entry
    def __enter__(self):
        # Save the default SIGCHLD handler to restore later
        self.default_SIGCHLD_handler = signal.getsignal(signal.SIGCHLD)
        return self

    # Executed upon context exit or exception
    def __exit__(self, exception_type, exception_value, exception_traceback):
        # Remove child processes starting with lowest descendents first
        for subproc in reversed(psutil.Process(os.getpid()).children(recursive=True)):
            subproc.terminate()

        # Restore original stdout/stderr system attributes
        sys.stdout.flush()
        sys.stderr.flush()
        for stream, orig_target in self._old_targets.items():
            setattr(sys, stream, orig_target)

        # Restore SIGCHLD handler to default
        signal.signal(signal.SIGCHLD, self.default_SIGCHLD_handler)

    def get_stdout(self):
        # If function_thread has been assigned an identity, then it has been .start()ed
        if self.function_thread.ident is not None:
            return self._stdout.getvalue()
        else:
            return ''

    def get_stderr(self):
        # If function_thread has been assigned an identity, then it has been .start()ed
        if self.function_thread.ident is not None:
            return self._stderr.getvalue()
        else:
            return ''

    # Check if the function has been executed AND is no longer alive
    def is_complete(self):
        return self.function_thread.ident is not None and not self.function_thread.is_alive()

    # Private method
    @staticmethod
    def __handle_sigchld(signum, stackframe):
        # Removes terminated child processes from the process table.
        # This avoids creating zombie processes after termination.
        try:
            os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            pass
