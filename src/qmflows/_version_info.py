from typing import NamedTuple, Final

from packaging.version import Version

from ._version import __version__

__all__ = ["version_info"]


class VersionInfo(NamedTuple):
    """A :func:`~collections.namedtuple` representing the version of a package."""

    #: :class:`int`: The semantic_ major version.
    major: int = 0

    #: :class:`int`: The semantic_ minor version.
    minor: int = 0

    #: :class:`int`: The semantic_ micro version.
    micro: int = 0

    @property
    def patch(self) -> int:
        """:class:`int`: An alias for :attr:`VersionInfo.micro`."""
        return self.micro

    @property
    def maintenance(self) -> int:
        """:class:`int`: An alias for :attr:`VersionInfo.micro`."""
        return self.micro

    @property
    def bug(self) -> int:
        """:class:`int`: An alias for :attr:`VersionInfo.micro`."""
        return self.micro


VERSION = Version(__version__)

version_info: Final = VersionInfo._make(VERSION.release[:3])
