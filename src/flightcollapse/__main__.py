"""Allow ``python -m flightcollapse`` as well as the console script.

Both forms should work, and until 0.2.0 only one did.  The difference matters
in a venv reached by absolute path rather than by activation -- which is how
every wrapper script around this package invokes it -- because
``$ENV/bin/python -m flightcollapse`` is the form that needs no PATH at all.
"""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
