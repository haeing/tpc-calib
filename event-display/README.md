# TPC Event Display

JSROOT-based event display for the E72 TPC analyzer data format.

The `include/` and `src/` directories are kept from `../../E72/Simul/tpc_analyzer_2025`.
The added converter writes ordinary ROOT objects that JSROOT can draw in a browser.
The 3D TPC frame is exported as `TPolyLine3D` wireframe objects, so it overlays
cleanly with hits and tracks in the same JSROOT view.
The browser does not load all entries at once. `make serve` starts a small local
server that converts only the requested entry when Previous/Next is pressed.

## Build

```sh
make
```

This requires ROOT and `root-config` in `PATH`.

## Create Display ROOT File

Geometry only:

```sh
make geometry
```

From an analyzer ROOT file:

```sh
make convert INPUT=data/evt_display.root TREE=tpc EVENT=0
```

Optional variables:

- `TREE=tree_name`
- `OUTPUT=web/display.root`
- `COUNT=1` or `COUNT=all` for offline batch export

## Open With JSROOT

```sh
make serve
```

Then open:

```text
http://localhost:8000/
```

In the web UI, set:

- input ROOT file path, for example `data/evt_display.root`
- tree name, for example `tpc` or `ktpc_g`
- entry number

Previous/Next calls the local render API and rewrites `web/display.root` with
only that entry.
