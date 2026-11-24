"""`python -m starnest run config.yaml` -> starnest.cli.main()."""
import sys

from .cli import main

sys.exit(main())
