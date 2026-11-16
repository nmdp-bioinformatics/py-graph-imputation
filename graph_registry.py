"""Thread-safe registry of named imputation graphs (name -> config file).

Shared by api_server.py and wmda_api_server.py. Graphs are built lazily on
first use, once per name, and reused across requests.
"""

import threading

from imputation_service import ImputationService


class GraphRegistry:
    """Registers named graphs (name -> config file) and lazily builds services."""

    def __init__(self):
        self._configs = {}
        self._services = {}
        self._default = None
        # Guards _configs/_services membership checks and the per-name lock map.
        self._registry_lock = threading.Lock()
        # One build lock per graph name, so concurrent first-time requests for
        # the SAME graph build it once instead of racing to build it N times.
        self._build_locks = {}

    def register(self, name, config_file, make_default=False):
        with self._registry_lock:
            self._configs[name] = config_file
            if make_default or self._default is None:
                self._default = name

    @property
    def names(self):
        return list(self._configs.keys())

    @property
    def default_name(self):
        return self._default

    def is_loaded(self, name):
        return name in self._services

    def config_file(self, name):
        return self._configs.get(name)

    def _build_lock_for(self, name):
        """Return (creating once) the dedicated build lock for a graph name."""
        with self._registry_lock:
            lock = self._build_locks.get(name)
            if lock is None:
                lock = threading.Lock()
                self._build_locks[name] = lock
            return lock

    def get(self, name=None):
        if name in (None, ""):
            name = self._default
        if name not in self._configs:
            raise KeyError(name)
        # Fast path: already built (dict reads are atomic under the GIL).
        service = self._services.get(name)
        if service is not None:
            return service
        # Slow path: serialize builds of THIS graph on its own lock so only the
        # first thread builds it; the rest wait and reuse the cached instance.
        with self._build_lock_for(name):
            service = self._services.get(name)
            if service is None:
                service = ImputationService(self._configs[name])
                self._services[name] = service
            return service

    def load_all(self):
        for name in self._configs:
            self.get(name)

    def discover(self, graphs_root="graphs"):
        """
        Register any completed graph found under graphs_root/<name>/config.json
        (a graph is "completed" when its csv/nodes.csv exists). Does not build
        the graph; registration is lazy. Returns the list of newly registered
        names.
        """
        from pathlib import Path
        import json as _json

        newly = []
        root = Path(graphs_root)
        if not root.is_dir():
            return newly
        for cfg in sorted(root.glob("*/config.json")):
            name = cfg.parent.name
            if name in self._configs:
                continue
            # Only register if the graph actually finished generating.
            try:
                conf = _json.loads(cfg.read_text())
                csv_dir = Path(conf.get("graph_files_path", cfg.parent / "csv"))
                nodes = csv_dir / conf.get("node_csv_file", "nodes.csv")
            except Exception:  # noqa: BLE001
                nodes = cfg.parent / "csv" / "nodes.csv"
            if Path(nodes).is_file():
                self.register(name, str(cfg))
                newly.append(name)
        return newly
