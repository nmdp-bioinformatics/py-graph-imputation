#!/bin/bash
#
# Multi-worker launch for the HLA Imputation API via gunicorn.
#
# Why gunicorn: imputation is CPU-bound (pure-Python phasing) and the GIL means
# threads don't parallelize it — you need multiple *processes*. gunicorn runs N
# worker processes so N /impute requests run truly in parallel.
#
# --preload + GRIM_EAGER: the graphs are registered and built ONCE in the
# gunicorn master process (at module import, via _bootstrap_from_env), then the
# workers are forked. Forked workers share that graph memory copy-on-write, so
# you don't pay Nx RAM or Nx the (slow) first-load, and no request eats the
# cold-start latency.
#
# Tunables (override via env):
#   WORKERS   number of worker processes (default 6; box has 10 cores)
#   PORT      bind port (default 5000)
#   TIMEOUT   worker timeout secs — must exceed the eager graph-build time
#             (large 4-locus graphs can take ~60-90s to load) (default 120)
#   LOGLEVEL  app log level (default INFO)

cd "$(dirname "$0")"  # run from the repo root regardless of the caller's cwd

WORKERS="${WORKERS:-6}"
PORT="${PORT:-5000}"
TIMEOUT="${TIMEOUT:-120}"
LOGLEVEL="${LOGLEVEL:-INFO}"

# Stop any running server (gunicorn masters/workers OR a plain api_server.py),
# regardless of how it was launched. The [x] trick avoids matching this script.
pkill -f "[g]unicorn.*api_server" 2>/dev/null
pkill -f "[p]ython3 api_server.py" 2>/dev/null

# Belt-and-suspenders: free the port if a stray process still holds it.
lsof -tiTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null | xargs -r kill 2>/dev/null

sleep 1

source venv/bin/activate

# Eager-preload graphs in the master before forking workers.
export GRIM_EAGER=1
export GRIM_GRAPHS_DIR="graphs"   # auto-discover finished graphs under graphs/
export GRIM_LOG_LEVEL="${LOGLEVEL}"

# Limit which graphs eager-load (the rest stay lazy — built on first request).
# This matters because each graph is large/slow to load; loading all of them at
# startup can be very slow and memory-heavy.
#   EAGER_LIMIT   eager-load at most N graphs (e.g. EAGER_LIMIT=2 ./script.sh)
#   EAGER_GRAPHS  comma-separated graph names to eager-load (most precise)
[ -n "${EAGER_LIMIT:-}" ]  && export GRIM_EAGER_LIMIT="${EAGER_LIMIT}"
[ -n "${EAGER_GRAPHS:-}" ] && export GRIM_EAGER_GRAPHS="${EAGER_GRAPHS}"

# To pin specific graphs instead of discovery, uncomment:
# export GRIM_GRAPHS="cau=conf/minimal-configuration.json,us=graphs/us_broad_4loc/config.json"
# export GRIM_NO_DISCOVER=1

exec gunicorn \
  --workers "${WORKERS}" \
  --preload \
  --timeout "${TIMEOUT}" \
  --bind "0.0.0.0:${PORT}" \
  --access-logfile - \
  --error-logfile - \
  'api_server:app'
