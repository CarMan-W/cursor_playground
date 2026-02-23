#!/usr/bin/env python3
"""Generate a watertight STL mesh for a solid sphere.

The mesh is centered at the origin and uses millimeters as units.
Default output is a 3 mm diameter sphere.
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    radius = args.diameter / 2.0
    triangles = build_sphere_mesh(
        radius_mm=radius,
        lat_segments=args.lat_segments,
        lon_segments=args.lon_segments,
    )
    write_ascii_stl(args.output, "sphere_3mm_solid", triangles)
    print(
        f"Generated '{args.output}' with diameter={args.diameter} mm, "
        f"facets={len(triangles)}."
    )


if __name__ == "__main__":
    main()
