#!/bin/bash
# Quick test script to verify the web API setup

set -e

cd "$(dirname "$0")"  # run from the repo root regardless of the caller's cwd

echo "🧪 Testing Web API Setup"
echo "======================="

# Activate venv
echo "📦 Activating virtual environment..."
source venv/bin/activate

# Test imports
echo "✓ Testing Python imports..."
python3 -c "import streamlit; print('  ✅ Streamlit')";
python3 -c "import flask; print('  ✅ Flask')";
python3 -c "import flask_cors; print('  ✅ Flask-CORS')";

# Check files exist
echo "✓ Checking files..."
test -f streamlit_config_ui.py && echo "  ✅ streamlit_config_ui.py"
test -f api_server.py && echo "  ✅ api_server.py"
test -f launch_web_api.sh && echo "  ✅ launch_web_api.sh"
test -f WEB_API_README.md && echo "  ✅ WEB_API_README.md"
test -f requirements.txt && echo "  ✅ requirements.txt updated"

# Check config
echo "✓ Checking configuration..."
test -f conf/minimal-configuration.json && echo "  ✅ conf/minimal-configuration.json"

# Check data
echo "✓ Checking data directory..."
test -d data/freqs && echo "  ✅ data/freqs directory exists"

echo ""
echo "✅ All checks passed!"
echo ""
echo "Next steps:"
echo "1. Launch the UI:        source venv/bin/activate && python3 streamlit_config_ui.py"
echo "2. Or use the launcher:  ./launch_web_api.sh ui"
echo "3. After config saved:   python3 api_server.py --config config-name.json"
echo ""
echo "📚 See WEB_API_README.md for full documentation"
