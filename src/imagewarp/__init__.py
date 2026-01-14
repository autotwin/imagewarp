# Placeholder for module discovery
import os
import matplotlib

# Allow manual override via MPLBACKEND environment variable
backend = os.environ.get("MPLBACKEND")
if backend:
    matplotlib.use(backend)
elif os.environ.get("DISPLAY"):
    matplotlib.use("TkAgg")  # Interactive backend for local machine
else:
    matplotlib.use("Agg")  # Non-interactive backend for server
