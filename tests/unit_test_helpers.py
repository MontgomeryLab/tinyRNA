"""
Helper functions for unit tests
"""

import os
import hashlib

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