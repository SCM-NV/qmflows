"""A module containing the :class:`~logging.Logger` of QMFlows."""

import sys
import logging

__all__ = ['logger', 'stdout_handler']

#: The QMFlows :class:`~logging.Logger`.
logger = logging.getLogger(__package__)
logger.setLevel(logging.DEBUG)

stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_handler.setLevel(logging.DEBUG)
stdout_handler.setFormatter(logging.Formatter(
    fmt='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
))
logger.addHandler(stdout_handler)
