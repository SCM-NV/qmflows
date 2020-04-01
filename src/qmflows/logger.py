"""A module containing the :class:`~logging.Logger` of QMFlows."""

import logging

__all__ = ['logger']

#: The QMFlows :class:`~logging.Logger`.
logger = logging.getLogger(__package__)
