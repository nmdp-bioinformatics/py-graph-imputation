# Web API Setup for HLA Imputation

This branch contains a web-based configuration UI and API server for the HLA imputation system.

## Quick Start

### 1. Install Dependencies

```bash
pip install streamlit flask flask-cors
```

### 2. Launch the Configuration UI

```bash
streamlit run streamlit_config_ui.py
```

This opens a web interface at `http://localhost:8501` where you can:
- Select populations from available frequency files
- Configure input/output loci
- Set algorithm parameters (epsilon, num_results, priority, plan B, etc.)
- Choose phased or unphased mode
- Save the configuration to `conf/`

### 3. Launch the API Server

After saving a configuration from the UI:

```bash
python api_server.py --config minimal-configuration.json --port 5000
```

Or with debug mode:

```bash
python api_server.py --config minimal-configuration.json --port 5000 --debug
```

## API Endpoints

- **GET `/`** - API documentation
- **GET `/health`** - Health check
- **GET `/config`** - Get current configuration
- **POST `/impute`** - Perform imputation
  - Request: `{"patient": "...", "donor": "..."}`
- **POST `/match`** - Legacy match endpoint (for compatibility)
  - Request: `{"patient_slug": "...", "donor_slug": "..."}`

## Workflow

1. **Configure**: Use Streamlit UI to set parameters
2. **Save**: Save configuration to `conf/` directory
3. **Launch**: Start API server with `api_server.py --config <config-name>.json`
4. **Use**: Send requests to API endpoints

## Configuration Parameters

The UI exposes the following parameters:

- **Populations**: Select from available frequency files in `data/freqs/`
- **Input/Output Loci**: Choose which loci to impute
- **Phasing**: Phased or unphased mode
- **Number of Results**: Max imputation results
- **Epsilon (ε)**: Convergence threshold
- **Frequency Trim Threshold**: Filter low frequencies
- **Plan B**: Enable fallback strategy
- **Priority Weights**: Alpha, beta, gamma, delta, eta

## Project Structure

```
streamlit_config_ui.py   # Web UI for configuration
api_server.py            # Flask API server
conf/                    # Configuration files (generated)
data/freqs/              # Population frequency files
data/subjects/           # Subject data
output/                  # Output files
```

## Notes

- Populations are auto-discovered from `data/freqs/`
- In the future, this will integrate with HFX/PHYCUS for remote frequency data
- Configuration files are stored as JSON in `conf/`
- The API server can be configured independently without the UI
