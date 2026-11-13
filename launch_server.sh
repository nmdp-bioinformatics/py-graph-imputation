#!/bin/bash

cd "$(dirname "$0")"  # run from the repo root regardless of the caller's cwd

# Kill any running api_server.py (matches however it was launched: with or
# without --graph/--config). -f matches the full command line; the pattern is
# anchored on the script name so it won't match this script or the editor.
pkill -f "[p]ython3 api_server.py" 2>/dev/null

# Belt-and-suspenders: free port 5000 in case a stray process still holds it.
lsof -tiTCP:5000 -sTCP:LISTEN 2>/dev/null | xargs -r kill 2>/dev/null

# Give the OS a moment to release the socket before rebinding.
sleep 1

source venv/bin/activate
python3 api_server.py --port 5000

# or log to a file
#nohup python3 api_server.py --graph cau=conf/minimal-configuration.json --port 5000 > /tmp/api_server.log 2>&1 &

exit 0;

# test
 curl -s -X POST http://localhost:5000/impute \
   -H 'Content-Type: application/json' \
   -d '{"graph":"us_broad_4loc","glstring":"A*01:XX^B*08:XX","output":"umug","marginal_ratio_threshold":0.001}' | python3 -m json.tool
