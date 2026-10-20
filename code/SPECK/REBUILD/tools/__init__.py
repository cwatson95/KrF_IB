"""REBUILD.tools — scripts allowed to import OTHER packages.

Never imported by MODULES/ (gated in tests/test_visualize_3d.py); the
coupling between packages is a tool's OUTPUT (npz, PNG, HTML), never an
import the engine sees.
"""
