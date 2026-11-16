# HLA Imputation API Server

A lightweight web service that performs single-input HLA imputation over a
frequency-backed graph, and can generate new graphs on demand. It wraps the
`grim` imputation engine so you can go from a GL string to ranked phased /
unphased genotypes and populations over HTTP, without the file-heavy config
workflow.

- Server entry point: `api_server.py`
- Imputation core (loads graph once, imputes in memory): `imputation_service.py`
- Graph generation factory: `graph_factory.py`
- Interactive API docs (Swagger UI): `http://<host>:<port>/apidocs`

---

## Quick start

```bash
# single graph, lazy-loaded on first request
python3 api_server.py --graph cau=conf/minimal-configuration.json --port 5000

# open the Swagger UI
open http://localhost:5000/apidocs
```

Impute (the new sample):

```bash
curl -X POST http://localhost:5000/impute \
  -H 'Content-Type: application/json' \
  -d '{
        "glstring": "A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01",
        "output": "both",
        "marginal_ratio_threshold": 0.001
      }'
```

---

## Concepts

### GL string format
- Loci separated by `^`
- The two haplotypes of a locus separated by `+`
- Allele ambiguity separated by `/`

Example: `A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01`

The engine itself only understands explicit allele-level ambiguity (`/`). It
does **not** natively expand ambiguity *codes* — but the server will, via
py-ard, before imputing. See
[Ambiguity codes (XX / MAC / serology)](#ambiguity-codes-xx--mac--serology).

### Graphs
A **graph** is a frequency-backed data structure built for a specific set of
populations and loci. Each graph is identified by a **name** and described by a
config JSON. The server:

- loads a graph **once** (lazily on first use, or eagerly with `--eager`) and
  reuses it across requests;
- **auto-detects** the graph's populations (from `info_node.csv`) and loci (from
  `nodes.csv`), so you don't have to keep the config's `populations` / `loci_map`
  in sync with the actual graph — the service adapts to whatever is loaded.

### Output flavors
- `umug` — **U**nphased **M**ulti-locus **U**nambiguous **G**enotype
  (internally the engine calls this MUUG).
- `pmug` — **P**hased **M**ulti-locus **U**nambiguous **G**enotype.
- `both` — return both sections.

Each result carries a raw `frequency` and a renormalized `probability`
(see below).

### Result truncation & renormalization (`marginal_ratio_threshold`)
Instead of (or in addition to) a fixed `number_of_results` cap, you can pass a
`marginal_ratio_threshold` (e.g. `0.001`). Results are ranked by frequency, then
walked from the top while accumulating a running sum. At the first result where

```
next_frequency / cumulative_sum_so_far  <  marginal_ratio_threshold
```

the list is truncated. The kept results are then **renormalized** so their
`probability` values sum to 1. This gives a natural, data-driven cutoff: the
list stops once additional results contribute negligibly to the total mass.

The threshold can be set per-request or as a `marginal_ratio_threshold` key in
the graph's config JSON.

### Ambiguity codes (XX / MAC / serology)

The imputation engine only understands explicit allele-level ambiguity via the
`/` separator (e.g. `A*01:01/A*01:02`). It has **no** built-in knowledge of
ambiguity *codes*:

- **XX codes** — `A*01:XX` (all alleles sharing the first field).
- **MAC codes** — `A*01:AB` (Multiple Allele Codes, e.g. NMDP codes).
- **Serology** — `A1`, `Bw4`, `Cw7`, `DR15`, …

To make these work, the server expands them with
[`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) **before** imputing,
using py-ard's `lgx` reduction. `lgx` both expands the code into its member
alleles and reduces every allele to the **2-field** ARD level (e.g.
`A*01:01`) — the same resolution the graphs are built at — producing a plain
`/`-ambiguous GL string the engine can impute over.

Controlled by the `expand` parameter:

- `auto` **(default)** — expand only if the input looks like it contains a
  code (a plain 2-field or `/`-ambiguous string is imputed directly, so it
  never pays the py-ard cost).
- `true` — always run py-ard's `lgx` redux on the input.
- `false` — never expand; pass the GL string through unchanged.

When an expansion happens, the response reflects it with a flag (the expanded
and original GL strings are **not** returned — an expanded XX code can be
enormous; the strings are written to the server log instead):

```jsonc
{
  "graph": "cau",
  "expanded": true,          // py-ard expanded a code in the input
  "status": "ok",
  "populations": ["CAU", "HIS", "NAM", "AFA", "API"],
  "umug": { ... }
}
```

#### Single-allele (homozygous) loci

The engine's genotype parser requires **two haplotypes per locus** (an
`allele+allele` genotype). A locus written as a single allele — e.g.
`A*01:01^B*08:01`, or a single code like `A*01:XX^B*08:XX` — is treated as
**homozygous**: the service duplicates the (possibly ambiguous) allele into both
haplotypes (`X` → `X+X`) before imputing. This normalization runs *after* code
expansion, so `A*01:XX^B*08:XX` becomes a homozygous, `/`-ambiguous genotype at
each locus. When it happens the response carries `"homozygous_filled": true`.
Loci that already have a `+` are left untouched.

Notes and caveats:

- **First init is slow.** py-ard's first `init()` downloads a large
  IPD-IMGT/HLA reference database (cached under `$TMPDIR/pyard-$USER/` by
  default). It is initialized **lazily** and cached, so only requests that
  actually contain codes trigger it.
- **Graceful degradation.** If py-ard is not installed, or it cannot expand the
  input (e.g. an invalid allele/code), the original string is imputed as-is and
  the response carries an `expand_note` explaining why nothing was expanded.
- **Pin the DB version** with an `imgt_version` key in the graph's config JSON
  if you need reproducible expansions across IPD-IMGT/HLA releases.

Example:

```bash
curl -X POST http://localhost:5000/impute \
  -H 'Content-Type: application/json' \
  -d '{
        "graph": "us_broad_4loc",
        "glstring": "A*01:XX+A*02:XX^B*07:XX+B*08:XX^C*07:XX+C*07:XX^DRB1*03:XX+DRB1*15:XX",
        "output": "both",
        "marginal_ratio_threshold": 0.001
      }'
```

---

## Endpoints

| Method | Path | Description |
|-------:|------|-------------|
| GET | `/` | Redirects to the Swagger UI (`/apidocs`). |
| GET | `/apidocs` | Swagger / OpenAPI UI. |
| GET | `/health` | Readiness plus which graphs are registered/loaded. |
| GET/POST | `/impute` | Impute a single GL string. |
| GET | `/graphs` | List registered graphs and their load state. |
| POST | `/graphs` | Generate + register a new graph (the factory). |
| GET | `/graphs/{name}` | Status of a graph (generation job + registration). |
| GET | `/factory/sources` | Frequency sources and the populations each offers. |
| GET | `/factory/suggest-name` | Suggest a graph name for a population/loci set. |

### `POST /impute`
Body (or query string) parameters:

| Field | Type | Notes |
|-------|------|-------|
| `glstring` | string | **Required.** Also accepts `gl_string` / `gl`. |
| `graph` | string | Which registered graph to use. Default = server default. |
| `pop1`, `pop2` | string | Optional population priors. |
| `output` | enum | `umug`, `pmug`, or `both` (default). |
| `number_of_results` | int | Cap on ranked genotype results. |
| `number_of_pop_results` | int | Cap on ranked population results. |
| `marginal_ratio_threshold` | float | Truncate + renormalize (see above). |
| `expand` | `auto`\|`true`\|`false` | Expand XX/MAC/serology codes via py-ard before imputing. Default `auto` (expand only when codes are detected). |

Response (abbreviated, `output=both`):

```json
{
  "graph": "cau",
  "populations": ["CAU", "HIS", "NAM", "AFA", "API"],
  "status": "ok",
  "umug": {
    "genotypes": [
      {"rank": 0, "genotype": "...", "frequency": 0.00116, "probability": 0.3468}
    ],
    "populations": [
      {"rank": 0, "pop": "CAU,CAU", "frequency": 0.9e-15, "probability": 1.0}
    ]
  },
  "pmug": {
    "genotypes": [
      {"rank": 0, "genotype": "hap1+hap2", "frequency": 0.00112, "probability": 0.3399}
    ],
    "populations": [ ... ]
  }
}
```

The response does **not** echo the input or imputed GL string. Additional
top-level fields appear only when relevant:

| Field | When present | Meaning |
|-------|--------------|---------|
| `status` | always | `ok` (results), `miss` (imputed to nothing), or `problem` (couldn't be processed). |
| `graph` | always | Name of the graph used. |
| `populations` | always | Populations auto-detected from the loaded graph. |
| `expanded` | code expansion happened | py-ard expanded an XX/MAC/serology code in the input. |
| `homozygous_filled` | single-allele locus filled | A locus written as one allele was treated as homozygous (`X` → `X+X`). |
| `expand_note` | expansion skipped/failed | Why the input was imputed unchanged (py-ard missing, invalid code, …). |
| `input_pops` | `pop1`/`pop2` supplied | Echoes the population priors used. |

---

## The graph factory

Generating a graph is a two-step pipeline (mirrors the main README):

1. `produce_hpf` — reads `<freq_data_dir>/<POP>.freqs.gz` and writes `hpf.csv`
   plus a population counts file.
2. `generate_graph` — reads `hpf.csv` and writes `nodes.csv`, `edges.csv`,
   `top_links.csv`, `info_node.csv`.

The factory runs both into an **isolated directory per graph** so nothing
collides:

```
graphs/<name>/
  config.json
  csv/
    hpf.csv
    pop_counts_file.txt
    nodes.csv, edges.csv, top_links.csv, info_node.csv
```

When generation finishes successfully, the graph is **auto-registered** and
immediately available to `/impute?graph=<name>`.

### Frequency sources / loci
The factory reads from **HFX bundles** (structured, self-describing `.hfx`
archives read via [`hfx-tools`](https://github.com/immunomath/hfx-tools)) or the
legacy `.freqs.gz` files. Each source is selected by a `source` id (or, as a
convenience, by `loci` count):

| `source` | Type | Loci | Dir | Populations |
|----------|------|------|-----|-------------|
| `hfx4` **(default)** | HFX | A, B, C, DRB1 | `~/src/HF/HFX/4loc` | AFA, CAU, HIS, NAM, API, … (26) |
| `freqs5` | freqs.gz | A, C, B, DRB1, DQB1 | `data/freqs` | CAU |
| `freqs6` | freqs.gz | A, C, B, DRB1, DQB1, DPB1 | `data/freqs/6loc` | CAU, HIS, NAM, AFA, API |

For HFX sources, the loci order and names are read directly from each bundle's
`metadata.json` (`outputResolution`), so the graph always matches the data.

Use `GET /factory/sources` to see exactly which populations are available in
each source on your machine.

### Create a graph

```bash
# 4-locus HFX CAU (default source). Name is optional — omit it and the
# server suggests one (e.g. "cau_hfx4").
curl -X POST http://localhost:5000/graphs \
  -H 'Content-Type: application/json' \
  -d '{"populations": ["CAU"], "source": "hfx4", "freq_trim_threshold": 0.00001}'
```

**Multi-population example** — a single graph over five populations
(the frequency vector on each node has one entry per population):

```bash
curl -X POST http://localhost:5000/graphs \
  -H 'Content-Type: application/json' \
  -d '{
        "name": "us5_hfx4",
        "populations": ["AFA","CAU","HIS","NAM","API"],
        "source": "hfx4",
        "freq_trim_threshold": 0.00001
      }'
```

Response (`202 Accepted`): generation runs in the background. Poll it:

```bash
curl http://localhost:5000/graphs/us5_hfx4
# {"name":"us5_hfx4","registered":false,"loaded":false,
#  "job":{"status":"running", ...}}   # -> "ready" when done, then registered:true
```

Once `registered:true`, impute against it with `"graph": "us5_hfx4"`.

`freq_trim_threshold` controls size/speed: **higher = fewer haplotypes = smaller
and faster** to build. Small single-population 4-locus HFX graphs build in
seconds; large populations or many-population graphs take longer.

### Suggested graph names
Naming graphs by hand is easy to get inconsistent, so the tool suggests names
for you. The convention is:

```
<populations joined by '-'>_<source suffix>
```

where the source suffix is the HFX source id (e.g. `hfx4`) or `<loci>loc` for
the legacy freqs sources. Many populations collapse to a count-based short
form, and a `-2`, `-3`, … suffix guarantees uniqueness.

| Populations | Source | Suggested name |
|-------------|--------|----------------|
| `["CAU"]` | `hfx4` | `cau_hfx4` |
| `["CAU","HIS"]` | `hfx4` | `cau-his_hfx4` |
| `["AFA","CAU","HIS","NAM","API"]` | `hfx4` | `afa-plus4_hfx4` |
| `["CAU"]` (legacy) | `freqs6` | `cau_6loc` |
| `["CAU"]` (when `cau_hfx4` exists) | `hfx4` | `cau_hfx4-2` |

Preview a name without creating anything:

```bash
curl "http://localhost:5000/factory/suggest-name?populations=CAU,HIS&source=hfx4"
# {"suggested_name":"cau-his_hfx4","populations":["CAU","HIS"],"source":"hfx4","loci":null}
```

If you `POST /graphs` without a `name`, the server applies this same suggestion
automatically and returns the chosen name in the job.

---

## Running options

```
python3 api_server.py [options]

--config PATH            Register a single config as graph "default".
--graph NAME=CONFIG      Register a named graph. Repeatable.
--eager                  Build all registered graphs at startup
                         (default: lazy — each builds on first /impute).
--eager-graphs NAMES     Comma-separated graph names to eager-load; the rest
                         stay lazy (built on first request).
--eager-limit N          Eager-load at most N graphs (the first N registered).
--graphs-dir DIR         Directory scanned for existing graphs to
                         auto-register (default: graphs).
--no-discover            Do not auto-register completed graphs found in
                         --graphs-dir.
--host HOST              Default 0.0.0.0
--port PORT              Default 5000
--debug                  Flask debug mode (reloader disabled to avoid
                         building graphs twice).
--log-level LEVEL        Logging level (default: INFO). One of DEBUG, INFO,
                         WARNING, ERROR, CRITICAL.
```

> **No config required.** `--config` and `--graph` are both optional. On
> startup the server registers any `--graph`/`--config` you pass, then
> **auto-discovers** finished graphs under `--graphs-dir` (so graphs built via
> the factory or CLI are picked up automatically). If nothing is registered at
> all, it falls back to registering `conf/minimal-configuration.json` as
> `default`.

Examples:

```bash
# Multiple graphs, lazy
python3 api_server.py \
  --graph cau=conf/minimal-configuration.json \
  --graph wmda=conf/WMDA-minimal-configuration.json \
  --port 5000

# Preload everything at startup
python3 api_server.py --graph cau=conf/minimal-configuration.json --eager
```

---

## Scaling: multiple workers (gunicorn)

`python3 api_server.py` uses Flask's **development server**, which effectively
serves one request at a time. Since imputation is **CPU-bound** (pure-Python
phasing) and Python has the GIL, throughput is limited to a single core and
concurrent `/impute` calls queue behind each other.

For real concurrency, run under **gunicorn** with multiple **worker processes**
(threads won't help CPU-bound work — you need separate processes):

```bash
GRIM_EAGER=1 gunicorn \
  --workers 4 \
  --preload \
  --timeout 300 \
  --bind 0.0.0.0:5000 \
  'api_server:app'
```

Or just use the provided script (defaults: 6 workers, preload, discovery):

```bash
./launch_server_gunicorn.sh
# override: WORKERS=8 PORT=5000 TIMEOUT=300 ./launch_server_gunicorn.sh
```

### Why `--preload` + `GRIM_EAGER` matters here

A graph is large (hundreds of MB of CSV → multi-GB in memory) and slow to load.
Without care, each of N workers would build its **own** copy — Nx the RAM and
Nx the load time. Instead:

- **`--preload`** makes gunicorn import `api_server` in the **master** process
  before forking workers.
- **`GRIM_EAGER=1`** triggers an import-time bootstrap (`_bootstrap_from_env`)
  that registers and **builds the graphs once in the master**.
- Workers are then **forked**, sharing the already-built graph memory
  copy-on-write. No worker rebuilds, and no request pays the cold-start cost.

> The per-graph build lock only protects against duplicate builds **within** a
> process. Across gunicorn workers it does not apply — which is exactly why
> eager-preloading in the master (before fork) is the right pattern.

### Bootstrap environment variables

Used by the import-time bootstrap (i.e. under gunicorn `--preload`):

| Variable | Meaning |
|----------|---------|
| `GRIM_EAGER` | Truthy → run bootstrap and eager-load graphs in the master. Required to preload. |
| `GRIM_GRAPHS` | Comma-separated `NAME=CONFIG` specs (like repeated `--graph`). |
| `GRIM_CONFIG` | Single config registered as graph `default`. |
| `GRIM_GRAPHS_DIR` | Directory to auto-discover finished graphs from (default `graphs`). |
| `GRIM_NO_DISCOVER` | Truthy → skip auto-discovery. |
| `GRIM_EAGER_GRAPHS` | Comma-separated graph **names** to eager-load; the rest stay lazy. |
| `GRIM_EAGER_LIMIT` | Eager-load at most N graphs (the first N registered). |
| `GRIM_LOG_LEVEL` | Logging level (default `INFO`). |

### Sizing workers

- Rule of thumb for CPU-bound work: **`workers ≈ physical cores`**. Leave a
  little headroom for the OS (e.g. 4–8 on a 10-core box).
- **RAM:** on Linux, forked workers share the graph copy-on-write, so memory
  stays close to a single graph's footprint. On macOS, copy-on-write is less
  effective, so budget roughly `workers × graph_size` and keep the worker count
  modest. With ~2 GB per graph and 32 GB RAM, 4–8 workers is comfortable.
- **`--timeout`** must exceed the eager build time (large 4-locus graphs can
  take ~60–90 s); the default 30 s would kill the master mid-build.

---

## Notes & limitations

- **First load is slow.** Building a graph (whether at startup or lazily on the
  first request) reads large node/edge CSVs into memory. Subsequent requests are
  fast because the graph is cached.
- **Concurrent first-hits build once.** Lazy loading is guarded by a per-graph
  lock: if several requests for the same not-yet-loaded graph arrive at once,
  the graph is built a single time and shared, rather than each request building
  its own copy.
- **Progress bars are silenced.** The underlying `grim` engine prints tqdm
  progress bars ("Vertices/Edges Creation") while building a graph. The server
  sets `TQDM_DISABLE=1` at startup so these don't spam the logs; set
  `TQDM_DISABLE=0` in the environment if you want them back. (The CLI tools are
  unaffected and still show bars.)
- **Generation is CPU-heavy.** The `generate_graph` step enumerates locus
  combinations; full-population 6-locus graphs can take minutes. Raise
  `freq_trim_threshold` for smaller/faster graphs while experimenting.
- **Populations/loci are auto-detected** from the loaded graph files, so a
  config whose `populations`/`loci_map` don't match its graph will still work —
  the graph wins.
- Dependencies for the server: `flask`, `flask-cors`, `flasgger`, `py-ard`
  (for XX/MAC/serology expansion), and `gunicorn` (multi-worker deployment) —
  see `requirements.txt`.

## Logging

Both the server (`api_server` logger) and the imputation core
(`imputation_service` logger) emit through Python's standard `logging`, at the
level set by `--log-level` (default `INFO`), formatted as:

```
2026-09-29 00:33:25,149 INFO imputation_service: graph ready from config='graphs/us_broad_4loc/config.json'
```

Logs go to stderr, so redirecting the process output captures them alongside
Flask's access log:

```bash
python3 api_server.py --port 5000 --log-level INFO > /tmp/api_server.log 2>&1 &
```

What gets logged:

| Level | Event |
|-------|-------|
| INFO | Graph build start / ready; imputation `miss` (with the offending GL string). |
| WARNING | `problem` results, `expand_note` conditions, bad requests (400), py-ard init failure. |
| ERROR/exception | Graph build failures; unexpected imputation errors (500), with traceback. |

Because the full (possibly huge, post-expansion) GL strings are no longer
returned in the response, they are recorded here instead — so a `problem`/`miss`
log line includes both the caller's original input and the imputed string for
debugging.
