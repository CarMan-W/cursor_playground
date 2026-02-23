#!/usr/bin/env python3
"""Generate watertight STL meshes for spheres in millimeters.

Supports:
- Solid sphere
- Sphere with centered internal void (sphere/cube)
- 2D X/Y array layout with configurable center spacing
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import List, Sequence, Tuple

Vec3 = Tuple[float, float, float]


@dataclass(frozen=True)
class Triangle:
    v1: Vec3
    v2: Vec3
    v3: Vec3


def subtract(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def normalize(v: Vec3) -> Vec3:
    length = math.sqrt(dot(v, v))
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def average(*points: Vec3) -> Vec3:
    count = float(len(points))
    return (
        sum(p[0] for p in points) / count,
        sum(p[1] for p in points) / count,
        sum(p[2] for p in points) / count,
    )


def orient_outward(v1: Vec3, v2: Vec3, v3: Vec3) -> Triangle:
    """Ensure triangle winding produces outward normal for origin-centered sphere."""
    n = cross(subtract(v2, v1), subtract(v3, v1))
    c = average(v1, v2, v3)
    if dot(n, c) < 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def orient_inward(v1: Vec3, v2: Vec3, v3: Vec3) -> Triangle:
    """Ensure triangle winding produces inward normal for origin-centered meshes."""
    n = cross(subtract(v2, v1), subtract(v3, v1))
    c = average(v1, v2, v3)
    if dot(n, c) > 0.0:
        return Triangle(v1=v1, v2=v3, v3=v2)
    return Triangle(v1=v1, v2=v2, v3=v3)


def build_sphere_mesh(radius_mm: float, lat_segments: int, lon_segments: int) -> List[Triangle]:
    if radius_mm <= 0.0:
        raise ValueError("Radius must be > 0.")
    if lat_segments < 2:
        raise ValueError("lat_segments must be >= 2.")
    if lon_segments < 3:
        raise ValueError("lon_segments must be >= 3.")

    top: Vec3 = (0.0, 0.0, radius_mm)
    bottom: Vec3 = (0.0, 0.0, -radius_mm)

    rings: List[List[Vec3]] = []
    for i in range(1, lat_segments):
        theta = math.pi * i / lat_segments
        z = radius_mm * math.cos(theta)
        ring_radius = radius_mm * math.sin(theta)
        ring: List[Vec3] = []
        for j in range(lon_segments):
            phi = 2.0 * math.pi * j / lon_segments
            x = ring_radius * math.cos(phi)
            y = ring_radius * math.sin(phi)
            ring.append((x, y, z))
        rings.append(ring)

    triangles: List[Triangle] = []
    first_ring = rings[0]
    last_ring = rings[-1]

    # Top cap.
    for j in range(lon_segments):
        v2 = first_ring[(j + 1) % lon_segments]
        v3 = first_ring[j]
        triangles.append(orient_outward(top, v2, v3))

    # Middle strips.
    for ring_idx in range(len(rings) - 1):
        ring_a = rings[ring_idx]
        ring_b = rings[ring_idx + 1]
        for j in range(lon_segments):
            v00 = ring_a[j]
            v01 = ring_a[(j + 1) % lon_segments]
            v10 = ring_b[j]
            v11 = ring_b[(j + 1) % lon_segments]
            triangles.append(orient_outward(v00, v01, v11))
            triangles.append(orient_outward(v00, v11, v10))

    # Bottom cap.
    for j in range(lon_segments):
        v2 = last_ring[j]
        v3 = last_ring[(j + 1) % lon_segments]
        triangles.append(orient_outward(bottom, v2, v3))

    return triangles


def build_cube_cavity_mesh(edge_mm: float) -> List[Triangle]:
    if edge_mm <= 0.0:
        raise ValueError("Cube void size must be > 0.")

    half = edge_mm / 2.0
    faces: Sequence[Tuple[Vec3, Vec3, Vec3, Vec3]] = (
        # +X
        ((half, -half, -half), (half, half, -half), (half, half, half), (half, -half, half)),
        # -X
        ((-half, -half, -half), (-half, -half, half), (-half, half, half), (-half, half, -half)),
        # +Y
        ((-half, half, -half), (-half, half, half), (half, half, half), (half, half, -half)),
        # -Y
        ((-half, -half, -half), (half, -half, -half), (half, -half, half), (-half, -half, half)),
        # +Z
        ((-half, -half, half), (half, -half, half), (half, half, half), (-half, half, half)),
        # -Z
        ((-half, -half, -half), (-half, half, -half), (half, half, -half), (half, -half, -half)),
    )

    triangles: List[Triangle] = []
    for a, b, c, d in faces:
        triangles.append(orient_inward(a, b, c))
        triangles.append(orient_inward(a, c, d))
    return triangles


def validate_internal_void(outer_radius_mm: float, void_shape: str, void_size_mm: float) -> None:
    if void_shape == "none":
        return

    if void_size_mm <= 0.0:
        raise ValueError("void_size must be > 0 when using an internal void.")

    if void_shape == "sphere":
        void_radius = void_size_mm / 2.0
        if void_radius >= outer_radius_mm:
            raise ValueError(
                "Sphere void diameter must be smaller than the outer sphere diameter."
            )
        return

    if void_shape == "cube":
        half_diagonal = (void_size_mm * math.sqrt(3.0)) / 2.0
        if half_diagonal >= outer_radius_mm:
            raise ValueError(
                "Cube void is too large; cube must fit fully inside the outer sphere."
            )
        return

    raise ValueError(f"Unsupported void_shape: {void_shape}")


def build_sphere_with_internal_void_mesh(
    radius_mm: float,
    lat_segments: int,
    lon_segments: int,
    void_shape: str,
    void_size_mm: float,
) -> List[Triangle]:
    triangles = build_sphere_mesh(
        radius_mm=radius_mm,
        lat_segments=lat_segments,
        lon_segments=lon_segments,
    )

    validate_internal_void(radius_mm, void_shape, void_size_mm)

    if void_shape == "sphere":
        cavity = build_sphere_mesh(
            radius_mm=void_size_mm / 2.0,
            lat_segments=lat_segments,
            lon_segments=lon_segments,
        )
        for tri in cavity:
            triangles.append(Triangle(v1=tri.v1, v2=tri.v3, v3=tri.v2))
    elif void_shape == "cube":
        triangles.extend(build_cube_cavity_mesh(void_size_mm))

    return triangles


def validate_array_config(array_x: int, array_y: int, array_spacing_mm: float) -> None:
    if array_x < 1:
        raise ValueError("array_x must be >= 1.")
    if array_y < 1:
        raise ValueError("array_y must be >= 1.")
    if (array_x > 1 or array_y > 1) and array_spacing_mm <= 0.0:
        raise ValueError("array_spacing must be > 0 when array_x or array_y is greater than 1.")


def offset_point(v: Vec3, dx: float, dy: float, dz: float) -> Vec3:
    return (v[0] + dx, v[1] + dy, v[2] + dz)


def offset_triangle(tri: Triangle, dx: float, dy: float, dz: float) -> Triangle:
    return Triangle(
        v1=offset_point(tri.v1, dx, dy, dz),
        v2=offset_point(tri.v2, dx, dy, dz),
        v3=offset_point(tri.v3, dx, dy, dz),
    )


def tile_mesh_2d(
    base_triangles: Sequence[Triangle],
    array_x: int,
    array_y: int,
    array_spacing_mm: float,
) -> List[Triangle]:
    validate_array_config(array_x, array_y, array_spacing_mm)

    if array_x == 1 and array_y == 1:
        return list(base_triangles)

    triangles: List[Triangle] = []
    x_center = (array_x - 1) / 2.0
    y_center = (array_y - 1) / 2.0

    for yi in range(array_y):
        dy = (yi - y_center) * array_spacing_mm
        for xi in range(array_x):
            dx = (xi - x_center) * array_spacing_mm
            for tri in base_triangles:
                triangles.append(offset_triangle(tri, dx, dy, 0.0))

    return triangles


def triangle_normal(tri: Triangle) -> Vec3:
    raw_normal = cross(subtract(tri.v2, tri.v1), subtract(tri.v3, tri.v1))
    return normalize(raw_normal)


def write_ascii_stl(path: str, solid_name: str, triangles: Sequence[Triangle]) -> None:
    with open(path, "w", encoding="ascii") as f:
        f.write(f"solid {solid_name}\n")
        for tri in triangles:
            n = triangle_normal(tri)
            f.write(f"  facet normal {n[0]:.8f} {n[1]:.8f} {n[2]:.8f}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {tri.v1[0]:.8f} {tri.v1[1]:.8f} {tri.v1[2]:.8f}\n")
            f.write(f"      vertex {tri.v2[0]:.8f} {tri.v2[1]:.8f} {tri.v2[2]:.8f}\n")
            f.write(f"      vertex {tri.v3[0]:.8f} {tri.v3[1]:.8f} {tri.v3[2]:.8f}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a watertight ASCII STL sphere mesh in millimeters."
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=3.0,
        help="Sphere diameter in mm (default: 3.0).",
    )
    parser.add_argument(
        "--lat-segments",
        type=int,
        default=48,
        help="Number of latitude segments from top to bottom (default: 48).",
    )
    parser.add_argument(
        "--lon-segments",
        type=int,
        default=96,
        help="Number of longitude segments around axis (default: 96).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="sphere_3mm_solid.stl",
        help="Output STL file path.",
    )
    parser.add_argument(
        "--void-shape",
        type=str,
        choices=("none", "sphere", "cube"),
        default="none",
        help="Internal void shape: none, sphere, or cube (default: none).",
    )
    parser.add_argument(
        "--void-size",
        type=float,
        default=1.0,
        help=(
            "Internal void size in mm. Sphere uses diameter; cube uses edge length "
            "(default: 1.0). Ignored when --void-shape=none."
        ),
    )
    parser.add_argument(
        "--array-x",
        type=int,
        default=1,
        help="Number of objects along X axis (default: 1).",
    )
    parser.add_argument(
        "--array-y",
        type=int,
        default=1,
        help="Number of objects along Y axis (default: 1).",
    )
    parser.add_argument(
        "--array-spacing",
        type=float,
        default=100.0,
        help=(
            "Center-to-center spacing in mm used for the 2D array "
            "(default: 100.0)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    radius = args.diameter / 2.0
    base_triangles = build_sphere_with_internal_void_mesh(
        radius_mm=radius,
        lat_segments=args.lat_segments,
        lon_segments=args.lon_segments,
        void_shape=args.void_shape,
        void_size_mm=args.void_size,
    )
    triangles = tile_mesh_2d(
        base_triangles=base_triangles,
        array_x=args.array_x,
        array_y=args.array_y,
        array_spacing_mm=args.array_spacing,
    )
    solid_name = "sphere_with_internal_void" if args.void_shape != "none" else "sphere_solid"
    write_ascii_stl(args.output, solid_name, triangles)
    if args.void_shape == "none":
        void_desc = "internal_void=none"
    else:
        void_desc = f"internal_void={args.void_shape}({args.void_size} mm)"
    if args.array_x == 1 and args.array_y == 1:
        array_desc = "array=1x1"
    else:
        array_desc = (
            f"array={args.array_x}x{args.array_y}, "
            f"spacing={args.array_spacing} mm"
        )
    print(
        f"Generated '{args.output}' with diameter={args.diameter} mm, "
        f"{void_desc}, {array_desc}, objects={args.array_x * args.array_y}, "
        f"facets={len(triangles)}."
    )


if __name__ == "__main__":
    main()
