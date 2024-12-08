#!/bin/bash

export LD_LIBRARY_PATH=/app/lib64:/app/lib:$LD_LIBRARY_PATH
exec /app/bin/FlatFab "$@"

