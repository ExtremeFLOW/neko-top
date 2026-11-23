"""Stand-in for the part of pkg_resources that flinter uses.

flinter 0.7.1, its latest release, imports pkg_resources to locate its default
rule files. setuptools removed pkg_resources in version 82, so this module
provides the one function flinter needs with importlib.resources instead, and
setuptools is not needed at all. Put this directory on PYTHONPATH when running
flint.
"""

import importlib.resources


def resource_filename(package, resource):
    """Return the path of a file inside an installed package."""
    return str(importlib.resources.files(package).joinpath(resource))
