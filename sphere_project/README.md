# cursor_playground
cursor_playground

## Generate a 3 mm solid sphere STL

This repository includes a Python script to generate a watertight STL mesh
for a solid sphere (surface representation) in millimeter units.

It also includes a browser-based generator page where users can input sphere
diameter and download a generated STL file directly.

### Run

```bash
cd sphere_project
python3 generate_sphere_stl.py
```

By default this creates:

- `sphere_3mm_solid.stl`
- Diameter: `3.0 mm` (radius `1.5 mm`)

### Optional parameters

```bash
cd sphere_project
python3 generate_sphere_stl.py \
  --diameter 3.0 \
  --lat-segments 48 \
  --lon-segments 96 \
  --output sphere_3mm_solid.stl
```

## Web UI: input diameter and download STL

Open this page in a browser:

- `sphere_project/index.html`

The page allows the user to:

1. Enter sphere diameter in millimeters.
2. Click **Generate and Download STL**.
3. Download the generated ASCII STL file.

For a local web server, you can run:

```bash
cd sphere_project
python3 -m http.server 8000
```

Then open `http://localhost:8000/index.html`.
