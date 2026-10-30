#!/bin/bash

# This script runs inside your Docker container.
# Since the container already has Python 3.11 and your libraries,
# we just run your worker script directly.
mkdir -p results logs
python3 htc_worker.py "$@"
