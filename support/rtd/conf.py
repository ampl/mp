# Sphinx configuration for readthedocs.

import os

master_doc = 'index'
html_theme = 'sphinxdoc'
templates_path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '_templates')]
