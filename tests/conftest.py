"""
Put scripts/ on the import path so the tests can import the aravalli_wa package and paths.py
without the repository having to be installed.

Author: Bhavik Harish Lodhia, Curtin University
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
