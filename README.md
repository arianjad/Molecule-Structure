# Molecular-Structure
Python code written by Arian to calculate molecular energy levles and shifts


NOTE: if the code import is note working, copy and paste this to the beginning of the first cell:

import os
import sys

def add_to_sys_path(target_dir_name="Source Code"):
    current_dir = os.path.abspath('')
    while True:
        if target_dir_name in os.listdir(current_dir):
            path = os.path.join(current_dir, target_dir_name)
            if path not in sys.path:
                sys.path.append(path)
                print(f"Added {path} to sys.path")
            break  # Exit the loop after adding the path
        parent_dir = os.path.dirname(current_dir)
        if current_dir == parent_dir:  # Reached the top of the directory tree
            print(f"Directory named '{target_dir_name}' not found.")
            break  # Exit the loop if no directory is found
        current_dir = parent_dir  # Move up one directory level

add_to_sys_path()  # Call the function to add path