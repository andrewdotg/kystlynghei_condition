"""aef_geobox.py - lattice-snapped GeoBox builder for AEF AOI queries.

Adapted from a fork of aef_loader by Geethen Singh (upstream:
jakenotjay/aef-loader, commit feff13e, Apache License 2.0). The published
`aef-loader` package does not include this function, so it is included here.
"""

from __future__ import annotations

import math

from odc.geo.geobox import GeoBox


def aoi_geobox(
    bbox: tuple[float, float, float, float],
    crs: str,
    resolution: float,
    bbox_crs: str | None = None,
    snap: bool = True,
) -> GeoBox:
    """Build a target ``GeoBox`` for an AOI, snapped to a global pixel lattice.

    ``GeoBox.from_bbox`` anchors the grid at the AOI's own corner, so two AOIs
    reprojected independently (e.g. adjacent tiles, or the same area at different
    times) land on *different* pixel grids and cannot be mosaicked without
    resampling. Snapping the origin to integer multiples of ``resolution``
    (anchored at 0, 0) makes every AOI at a given resolution/CRS share one grid,
    so outputs align exactly and can be merged losslessly.

    Args:
        bbox: AOI bounds ``(minx, miny, maxx, maxy)``.
        crs: Target CRS for the GeoBox (e.g. ``"EPSG:32633"``).
        resolution: Pixel size in target CRS units (metres for UTM).
        bbox_crs: CRS of ``bbox`` if different from ``crs``; the bbox is
            reprojected (densified) to ``crs`` first. Defaults to ``crs``.
        snap: When True (default), snap the origin to the ``resolution`` lattice
            and grow the extent outward to fully cover ``bbox``. When False,
            behaves like ``GeoBox.from_bbox`` (AOI-anchored grid).

    Returns:
        A north-up ``GeoBox`` in ``crs`` at ``resolution`` covering ``bbox``.
    """
    minx, miny, maxx, maxy = bbox
    if bbox_crs is not None and str(bbox_crs) != str(crs):
        from pyproj import Transformer

        transformer = Transformer.from_crs(bbox_crs, crs, always_xy=True)
        minx, miny, maxx, maxy = transformer.transform_bounds(
            minx, miny, maxx, maxy, densify_pts=21
        )

    if not snap:
        return GeoBox.from_bbox(
            (minx, miny, maxx, maxy), crs=crs, resolution=resolution
        )

    # Snap outward to the global lattice: floor the min edges, ceil the max edges.
    snapped_minx = math.floor(minx / resolution) * resolution
    snapped_miny = math.floor(miny / resolution) * resolution
    snapped_maxx = math.ceil(maxx / resolution) * resolution
    snapped_maxy = math.ceil(maxy / resolution) * resolution
    return GeoBox.from_bbox(
        (snapped_minx, snapped_miny, snapped_maxx, snapped_maxy),
        crs=crs,
        resolution=resolution,
    )
