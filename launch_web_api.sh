#!/bin/bash
# Launch script for HLA Imputation Web API

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}═══════════════════════════════════════════${NC}"
    echo -e "${BLUE}  HLA Imputation Web API Launcher${NC}"
    echo -e "${BLUE}═══════════════════════════════════════════${NC}"
}

print_step() {
    echo -e "\n${GREEN}→${NC} $1"
}

print_info() {
    echo -e "${YELLOW}ℹ${NC} $1"
}

show_usage() {
    cat << EOF
Usage: ./launch_web_api.sh [COMMAND]

Commands:
    ui              Launch configuration UI only (default)
    server          Launch API server (requires config)
    both            Launch both UI and server
    help            Show this help message

Examples:
    ./launch_web_api.sh ui
    ./launch_web_api.sh server
    ./launch_web_api.sh both

Environment Variables:
    CONFIG_FILE     Config file to use with server (default: minimal-configuration.json)
    PORT            Port for API server (default: 5000)
    DEBUG           Enable debug mode (default: false)

EOF
}

# Main script
print_header

if [[ $# -eq 0 ]]; then
    COMMAND="ui"
else
    COMMAND="$1"
fi

case "$COMMAND" in
    ui)
        print_step "Launching Configuration UI"
        print_info "Opening http://localhost:8501"
        python3 -m streamlit run streamlit_config_ui.py
        ;;

    server)
        CONFIG="${CONFIG_FILE:-minimal-configuration.json}"
        PORT="${PORT:-5000}"
        print_step "Launching API Server"
        print_info "Config: $CONFIG"
        print_info "Port: $PORT"
        python3 api_server.py --config "$CONFIG" --port "$PORT"
        ;;

    both)
        CONFIG="${CONFIG_FILE:-minimal-configuration.json}"
        PORT="${PORT:-5000}"
        print_step "Launching Configuration UI"
        print_info "Opening http://localhost:8501"
        print_step "After saving config, you can start the API server with:"
        print_info "python3 api_server.py --config $CONFIG --port $PORT"
        python3 -m streamlit run streamlit_config_ui.py
        ;;

    help|--help|-h)
        show_usage
        ;;

    *)
        echo "Unknown command: $COMMAND"
        show_usage
        exit 1
        ;;
esac
