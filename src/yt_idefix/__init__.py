# import the api so that yt.load can be used
# immediately after `import yt.extensions.idefix`
from yt_idefix.api import *


from importlib.metadata import version

__version__ = version("yt-idefix")

del version
