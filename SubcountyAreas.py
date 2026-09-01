"""
Subcounty Areas

Requires: geopandas, shapely, pyproj, requests, scipy, numpy, pandas
    pip install geopandas shapely pyproj requests scipy numpy pandas

Groups census tracts into "subcounty areas" - a geography coarser than a tract, finer
than a county - by repeatedly merging the CURRENT LOWEST-POPULATION area into its
most-similar touching neighbor - similarity judged on both population and
residential-structure pattern (see step 2) - until reaching each target count in
AGGREGATION_TARGET_COUNTS. Runs twice over: constrained to stay inside one county, and
unconstrained. Uses a single tract source and a single residential-structure-location
source, each configured below via a *_SERVICE_URL constant (see the CONFIGURATION
section for the specific sources this template ships with). Adapting this to a
different tract source or region means editing those two URLs, their matching
ID/population/county field names, PROJECTED_CRS_EPSG, and likely
NEIGHBORHOOD_DISTANCES_METERS (see each constant's comment below) - the merge logic
itself is not tied to any particular source or region. Output is a single GeoPackage
(.gpkg) on the current user's Desktop, plus a CSV logging every merge decision for
later review.

1. Download and filter the input data.
   What: Pull every residential structure location and the configured tract
   boundaries over the ArcGIS REST API, reprojected to a single projected coordinate
   system in meters (PROJECTED_CRS_EPSG); keep only residential sites, and treat a
   residential unit count of 0 as 1.
   Why: Residential structures are the proxy for where people actually live; a unit count
   of 0 almost always means missing data rather than zero housing, so treating it as 1 is
   the conservative assumption.

2. Compute a residential-structure pattern signal (max_diff_k) for every original tract.
   What: Run Ripley's K across six distances (500-3,000m by default, chosen as a practical
   "neighborhood" scale for this template's default study area - reconsider
   NEIGHBORHOOD_DISTANCES_METERS for an area with a different typical settlement scale)
   and keep whichever single value deviates furthest from what pure
   randomness would predict, sign included. Edge correction is done by "plus sampling":
   every residential site actually inside the area is used as a center (none are discarded
   near the boundary), and neighbor counts for each center are drawn from a wider search
   area - the area's own boundary buffered outward by the largest tested distance (3,000m)
   - so a center near the edge still sees its true neighbors just across the line.
   Why: Raw population density mostly just reflects whether a tract was originally drawn as
   urban or rural (tracts are sized to equalize population, not area, so rural ones are
   large and urban ones are small by construction); Ripley's K instead measures whether
   homes are bunched together or spread apart relative to random chance, independent of the
   tract's raw size. Areas with fewer than 3 intersecting residential sites get max_diff_k =
   0, since the statistic can't be computed on that few points.

3. Repeatedly merge the current lowest-population area into its most-similar touching
   neighbor.
   What: Each round, find every pair of areas that share a real border and rank every pair
   by how close their populations are and how close their max_diff_k values are: each of
   the two differences is converted to its RANK among every pair touching map-wide that
   round (1 = most similar on that measure), and the two ranks are summed - the same
   ranking used for any pair anywhere on the map. Then take the single current area with the
   lowest population as this round's merge "seed" (skipping, for that round only, any area
   with no viable touching neighbor - there is no persistent exclusion list, so a skipped
   area is reconsidered fresh every round), and merge it with whichever of its OWN touching
   neighbors ranks best on that same metric.
   Why: If merges were instead chosen by picking whichever touching pair anywhere on the
   map is the closest match each round, nothing would ever force attention to a specific
   low-population area - an outlier whose only neighbors are all dissimilar to it could be
   passed over indefinitely while better-matched pairs elsewhere keep getting merged first,
   leaving it stranded near the population floor for most of the run. Seeding on the lowest
   population every round instead guarantees that area's population problem is addressed
   every single round. The cost is that a seed is occasionally forced to merge with a
   touching neighbor that is a poor residential-structure pattern match, when that's the only neighbor
   available; this cost is measured, not assumed, via the match-quality percentile logged
   for every merge in step 6. Ranking each difference (rather than dividing it by the sum
   of all differences map-wide) is what keeps that measurement meaningful pair-by-pair: a
   single pair with an extreme population or clustering difference elsewhere on the map
   would otherwise compress every other pair's normalized share toward zero, drowning out a
   real difference between two pairs actually being compared.

4. Track which original tracts end up in which merged area.
   What: Maintain a crosswalk table (original tract -> current area ID), updated every time
   two areas merge, and save it as a plain (non-spatial) table in the output GeoPackage
   alongside every snapshot.
   Why: Nothing about an original tract's own data is thrown away or approximated except
   max_diff_k, so any other tract-level statistic can always be correctly re-aggregated
   later from the source data instead of trusting a shortcut computed during merging.

5. Snapshot the result at each target aggregation count - or the closest count actually
   reached, if a run can't get all the way down to a target - plus the original, unmerged
   tracts themselves as a "BaseTracts" snapshot.
   What: Whenever the area count matches one of AGGREGATION_TARGET_COUNTS (60, 89, 100,
   107, 150), save both the current area boundaries and the crosswalk as of that point. Any
   run can run out of eligible touching pairs before reaching every target - nothing
   guarantees the remaining areas stay mutually reachable once merges are restricted to
   stay within a county, which makes this far more likely for the county-constrained run,
   but the unconstrained run isn't structurally immune either (e.g. tracts on an island with
   no land border to the mainland); when that happens, the run's actual stopping point is
   saved anyway as the closest available stand-in for whichever target(s) it never reached,
   rather than silently producing nothing for that cutoff. Separately, before any merging
   starts, the untouched original tracts are saved as their own "BaseTracts" layer, in the
   same schema as every other snapshot.
   Why: Produces a family of alternative aggregation levels to choose from, rather than a
   single fixed answer - and makes sure the county constraint's cost is a
   smaller-than-hoped-for area count, never a missing one. The original tract count (well
   above the largest target) would otherwise never appear in the output at all, even though
   it's the natural finest-resolution comparison point for anything joined against this
   output later.

6. Log every merge decision for later review.
   What: For every merge, record the seed and its chosen partner, their population and
   max_diff_k difference, the combined similarity score, and - critically - where that
   chosen match ranks (by percentile) against every pair available anywhere on the
   map that same round (including the chosen pair itself).
   Why: "Poor match" is only meaningful relative to what else was available that round;
   logging the percentile directly, rather than judging against a fixed threshold, is what
   lets that trade-off be measured empirically instead of just theorized about, and lets the
   same measurement be repeated on any future run.
"""

from __future__ import annotations

import math
import os
import sqlite3
from dataclasses import dataclass

import numpy as np
import pandas as pd
import geopandas as gpd
import requests
import shapely
from scipy.spatial import distance as spatial_distance
from scipy.stats import rankdata

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

OUTPUT_GEOPACKAGE_PATH = os.path.join(
    os.path.expanduser("~"), "Desktop", "SubcountyAreas", "SubcountyAreas.gpkg"
)
MERGE_DIAGNOSTICS_CSV_PATH = os.path.join(
    os.path.expanduser("~"), "Desktop", "SubcountyAreas", "SubcountyAreas_MergeDiagnostics.csv"
)

# Every residential structure location statewide. Defaults to Vermont's E911
# site/structure layer (VCGI) - replace with the equivalent address-point/structure
# layer for a different state.
RESIDENTIAL_SITES_SERVICE_URL = (
    "https://services1.arcgis.com/BkFxaEFNwHqX3tAw/arcgis/rest/services/"
    "FS_VCGI_OPENDATA_Emergency_SiteStructureAddressPoint_point_SP_v1_VIEW/FeatureServer/0"
)
# The field in that source that identifies the structure's category (residential, commercial, etc.)
RESIDENTIAL_CATEGORY_FIELD = "Category"
RESIDENTIAL_UNIT_COUNT_FIELD = "ResUnitCnt"
RESIDENTIAL_CATEGORY_VALUE = "Residential"

# The single tract source this run aggregates, and the field names within it that
# identify each tract's ID, total population, and county. Defaults to Vermont's 2020
# census tract layer (VCGI). To point this at a different tract vintage or a
# different state's tracts entirely, replace TRACT_SERVICE_URL and all three field
# names below to match whatever that source actually provides.

# 2020 tracts
TRACT_SERVICE_URL = (
    "https://services1.arcgis.com/BkFxaEFNwHqX3tAw/arcgis/rest/services/"
    "FS_VCGI_Census2020Tract_WM_v1/FeatureServer/0"
)
# tract field from 2020 census layer
TRACT_POPULATION_FIELD = "P0010001"

# 2010 tracts
#TRACT_SERVICE_URL = (
#    "https://services1.arcgis.com/BkFxaEFNwHqX3tAw/arcgis/rest/services/"
#    "FS_VCGI_OPENDATA_Demo_TRACT2010_poly_SP_v1/FeatureServer/0"
#)
# tract field from 2010 census layer
#TRACT_POPULATION_FIELD = "P0030001"

TRACT_ID_FIELD = "TRACT"
COUNTY_FIELD = "COUNTY"

# "Counties" keeps every merge inside one county; "NoCounties" allows any touching
# pair - so a NoCounties-mode area's county output field can end up listing more
# than one county, comma-separated (see combine_county_values).
COUNTY_CONSTRAINED_MODE = "Counties"
COUNTY_UNCONSTRAINED_MODE = "NoCounties"
COUNTY_MODES = [COUNTY_CONSTRAINED_MODE, COUNTY_UNCONSTRAINED_MODE]

# Cluster counts at which the current aggregation state gets snapshotted as output.
# Chosen with Vermont's own scale in mind (roughly 60-150 areas statewide) -
# reconsider these for a state with a very different number of tracts.
AGGREGATION_TARGET_COUNTS = [60, 89, 100, 107, 150]

# A single projected CRS, in meters - every download is reprojected to this so that
# every distance, area, and perimeter calculation below is in real-world meters.
# Defaults to Vermont State Plane; replace with an appropriate projected CRS (in
# meters) for a different state or region.
PROJECTED_CRS_EPSG = 32145

# The practical "neighborhood" scale Ripley's K is scanned across - chosen for
# Vermont's settlement density. Reconsider these distances for an area with a very
# different typical scale of clustering (e.g. a dense urban region).
NEIGHBORHOOD_DISTANCES_METERS = [500, 1000, 1500, 2000, 2500, 3000]


# -----------------------------------------------------------------------------
# RECORD TYPES
# -----------------------------------------------------------------------------

@dataclass
class SubcountyArea:
    """One current area on the map - either an original tract, or a merged
    cluster of tracts, at some point during the hierarchical merge."""
    cluster_id: str
    shape: shapely.Geometry
    population: int
    county: str  # comma-separated, deduplicated county list once a merge spans more than one
    max_diff_k: float


@dataclass
class MergeCandidate:
    """A touching pair of areas, and how similar they are to one another."""
    area_id_a: str
    area_id_b: str
    population_difference: float
    clustering_difference: float
    similarity_score: float = 0.0


@dataclass
class ResidentialSites:
    """Every residential structure location, ready for spatial queries: an
    in-memory spatial index, and the same points' coordinates/unit-count
    weights in the matching array order the index reports positions in.
    These three always travel together, so they're bundled here rather than
    passed around as three separate parameters."""
    spatial_index: shapely.STRtree
    coordinates: np.ndarray
    weights: np.ndarray


# -----------------------------------------------------------------------------
# DATA ACCESS
# -----------------------------------------------------------------------------

def fetch_feature_server_as_geodataframe(
    service_url: str, where_clause: str = "1=1", page_size: int = 1000
) -> gpd.GeoDataFrame:
    """Download every feature from an ArcGIS FeatureServer/MapServer layer's REST
    query endpoint, paging through resultOffset until a page comes back short,
    and return everything as one GeoDataFrame in the configured projected CRS."""
    downloaded_pages = []
    records_downloaded = 0
    while True:
        response = requests.get(
            f"{service_url}/query",
            params={
                "where": where_clause,
                "outFields": "*",
                "f": "geojson",
                "resultOffset": records_downloaded,
                "resultRecordCount": page_size,
                "returnGeometry": "true",
                "outSR": PROJECTED_CRS_EPSG,
            },
            timeout=120,
        )
        response.raise_for_status()
        features = response.json().get("features", [])
        if not features:
            break
        downloaded_pages.append(gpd.GeoDataFrame.from_features(features))
        records_downloaded += len(features)
        if len(features) < page_size:
            break
    combined_features = pd.concat(downloaded_pages, ignore_index=True)
    return gpd.GeoDataFrame(combined_features, geometry="geometry").set_crs(
        epsg=PROJECTED_CRS_EPSG, allow_override=True
    )


def prepare_datasets() -> tuple[ResidentialSites, gpd.GeoDataFrame]:
    """Download the filtered residential site points and the configured tract
    source. Returns (residential_sites, tract_geodata)."""
    print("Downloading residential site points...")
    residential_site_geodata = fetch_feature_server_as_geodataframe(
        RESIDENTIAL_SITES_SERVICE_URL, where_clause=f"{RESIDENTIAL_CATEGORY_FIELD} = '{RESIDENTIAL_CATEGORY_VALUE}'"
    )
    # A unit count of zero is treated conservatively as at least one unit.
    residential_site_geodata[RESIDENTIAL_UNIT_COUNT_FIELD] = (
        residential_site_geodata[RESIDENTIAL_UNIT_COUNT_FIELD].fillna(0).replace(0, 1)
    )
    print(f"  {len(residential_site_geodata)} residential sites downloaded.")

    residential_sites = ResidentialSites(
        spatial_index=shapely.STRtree(residential_site_geodata.geometry.values),
        coordinates=np.array([[point.x, point.y] for point in residential_site_geodata.geometry]),
        weights=residential_site_geodata[RESIDENTIAL_UNIT_COUNT_FIELD].to_numpy(dtype=float),
    )

    print("Downloading tract boundaries...")
    tract_geodata = fetch_feature_server_as_geodataframe(TRACT_SERVICE_URL)
    print(f"  {len(tract_geodata)} tracts downloaded.")

    return residential_sites, tract_geodata


def points_within_polygon(
    residential_sites: ResidentialSites, polygon: shapely.Geometry,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (global_indices, coordinates, weights) for every residential site
    intersecting the given polygon. The global indices let callers match the same
    point across two different polygon queries (e.g. an area and its buffered
    neighborhood), which a plain coordinate/weight pair can't do reliably."""
    matching_indices = residential_sites.spatial_index.query(polygon, predicate="intersects")
    return (
        matching_indices,
        residential_sites.coordinates[matching_indices],
        residential_sites.weights[matching_indices],
    )


# -----------------------------------------------------------------------------
# RESIDENTIAL-STRUCTURE PATTERN (Ripley's K) MEASURE
# -----------------------------------------------------------------------------

def compute_max_diff_k(residential_sites: ResidentialSites, polygon: shapely.Geometry) -> float:
    """Weighted Ripley's K, collapsed to a single value: whichever of the six
    tested distances deviates furthest from what pure randomness would predict,
    sign included (positive = clustered, negative = dispersed). Every
    residential site intersecting the polygon is used as a center - none are
    discarded near the boundary - because neighbor counts are drawn from a
    wider search area (the polygon buffered outward by the largest tested
    distance), so a center's true neighborhood is backed by real data even when
    it's close to the polygon's edge (plus sampling). The one place this can't
    help is the outer edge of the downloaded residential-site data itself (e.g. a
    state boundary), where there's no residential site data beyond the line to
    draw on.

    The "what would pure randomness predict" baseline intensity is estimated
    from that same wider buffered region, not from the polygon's own area/count.
    Using the polygon's own density would implicitly assume it's representative
    out to the full tested distance, which fails specifically for small polygons
    (they're small precisely because they're locally denser than their
    surroundings) - it would read every small, dense tract as extremely
    "dispersed" at the larger distance bands, since no real surrounding area
    sustains that tract's own hyper-local density that far out."""
    own_indices, own_coordinates, own_weights = points_within_polygon(residential_sites, polygon)
    point_count = len(own_coordinates)
    if point_count < 3:
        return 0.0

    own_weight_total = own_weights.sum()

    neighborhood_polygon = polygon.buffer(max(NEIGHBORHOOD_DISTANCES_METERS))
    neighbor_indices, neighbor_coordinates, neighbor_weights = points_within_polygon(
        residential_sites, neighborhood_polygon
    )
    neighbor_position_by_global_index = {
        global_index: position for position, global_index in enumerate(neighbor_indices)
    }
    baseline_intensity = neighbor_weights.sum() / neighborhood_polygon.area

    pairwise_distances = spatial_distance.cdist(own_coordinates, neighbor_coordinates)
    for row_index, global_index in enumerate(own_indices):
        pairwise_distances[row_index, neighbor_position_by_global_index[global_index]] = np.inf  # a point is not its own neighbor

    diff_k_values = []
    for distance_band in NEIGHBORHOOD_DISTANCES_METERS:
        within_distance_band = pairwise_distances <= distance_band
        weighted_neighbor_sums = (within_distance_band * neighbor_weights).sum(axis=1)
        total_weighted_pairs = (weighted_neighbor_sums * own_weights).sum()

        ripleys_k = total_weighted_pairs / (own_weight_total * baseline_intensity)
        ripleys_l = math.sqrt(ripleys_k / math.pi) if ripleys_k > 0 else 0.0
        diff_k_values.append(ripleys_l - distance_band)

    strongest_clustering, strongest_dispersion = max(diff_k_values), min(diff_k_values)
    return strongest_clustering if abs(strongest_clustering) >= abs(strongest_dispersion) else strongest_dispersion


# -----------------------------------------------------------------------------
# HIERARCHICAL MERGE - lowest-population seed merges with its best touching
# match, ranked map-wide, with per-merge diagnostics
# -----------------------------------------------------------------------------

def find_touching_area_pairs(
    subcounty_areas: dict[str, SubcountyArea], county_mode: str
) -> list[MergeCandidate]:
    """Every pair of currently-touching subcounty areas, paired with how
    different they are in population and in clustering signature. Uses an
    in-memory spatial index (STRtree), so this is a direct geometry query
    rather than a repeated call out to an external geoprocessing tool."""
    area_ids = list(subcounty_areas.keys())
    area_geometries = [subcounty_areas[area_id].shape for area_id in area_ids]
    area_spatial_index = shapely.STRtree(area_geometries)

    candidate_pairs = []
    pairs_already_found = set()
    for area_position, area_id in enumerate(area_ids):
        touching_positions = area_spatial_index.query(area_geometries[area_position], predicate="touches")
        for neighbor_position in touching_positions:
            neighbor_area_id = area_ids[neighbor_position]
            if neighbor_area_id == area_id:
                continue
            if county_mode == COUNTY_CONSTRAINED_MODE and (
                subcounty_areas[area_id].county != subcounty_areas[neighbor_area_id].county
            ):
                continue

            pair_key = tuple(sorted((area_id, neighbor_area_id)))
            if pair_key in pairs_already_found:
                continue
            pairs_already_found.add(pair_key)

            area_a, area_b = subcounty_areas[pair_key[0]], subcounty_areas[pair_key[1]]
            candidate_pairs.append(
                MergeCandidate(
                    area_id_a=pair_key[0],
                    area_id_b=pair_key[1],
                    population_difference=abs(area_a.population - area_b.population),
                    clustering_difference=abs(area_a.max_diff_k - area_b.max_diff_k),
                )
            )
    return candidate_pairs


def rank_candidates_by_similarity(candidate_pairs: list[MergeCandidate]) -> list[MergeCandidate]:
    """Score every touching pair by how similar they are, and sort most to
    least similar. Each of the two difference measures is converted to its
    RANK among every pair touching map-wide this round (1 = most similar on
    that measure), and the two ranks are summed. Ranking only depends on
    relative order, so a single pair with an extreme population or clustering
    difference elsewhere on the map can't compress every other pair's score
    toward zero and erase a real difference between two pairs actually being
    compared - which normalizing by the map-wide total would do. This ranking
    stays map-wide (not seed-relative) even though only the seed's own
    neighbors are eligible to be chosen below - so a seed's best available
    match is ranked on the same scale as everything else on the map that
    round, which is what makes the diagnostic percentiles in
    run_hierarchical_merge meaningful."""
    population_ranks = rankdata([pair.population_difference for pair in candidate_pairs])
    clustering_ranks = rankdata([pair.clustering_difference for pair in candidate_pairs])
    for pair, population_rank, clustering_rank in zip(candidate_pairs, population_ranks, clustering_ranks):
        pair.similarity_score = float(population_rank + clustering_rank)
    return sorted(candidate_pairs, key=lambda pair: pair.similarity_score)


def select_lowest_population_seed(
    subcounty_areas: dict[str, SubcountyArea], areas_with_candidate_pairs: set[str]
) -> str:
    """The current lowest-population area, considering only areas that have at
    least one viable touching neighbor this round. An isolated area (no
    same-mode neighbor) is skipped for THIS round only - there is no
    persistent exclusion list, so it is reconsidered fresh every round."""
    viable_area_ids = [area_id for area_id in subcounty_areas if area_id in areas_with_candidate_pairs]
    return min(viable_area_ids, key=lambda area_id: subcounty_areas[area_id].population)


def build_initial_areas_and_crosswalk(
    tract_geodata: gpd.GeoDataFrame,
) -> tuple[dict[str, SubcountyArea], pd.DataFrame]:
    """One SubcountyArea per original tract, and a matching identity crosswalk
    (every tract mapped to itself) - the starting point every hierarchical merge
    run is built from, and also what gets saved directly as the "BaseTracts"
    snapshot before any merging happens."""
    tract_crosswalk = pd.DataFrame({
        "OriginalTract": tract_geodata[TRACT_ID_FIELD].astype(str),
        "ClusterID": tract_geodata[TRACT_ID_FIELD].astype(str),
    })

    subcounty_areas: dict[str, SubcountyArea] = {
        str(tract_row[TRACT_ID_FIELD]): SubcountyArea(
            cluster_id=str(tract_row[TRACT_ID_FIELD]),
            shape=tract_row.geometry,
            population=int(tract_row[TRACT_POPULATION_FIELD])
            if pd.notna(tract_row[TRACT_POPULATION_FIELD])
            else 0,
            county=tract_row[COUNTY_FIELD],
            max_diff_k=tract_row["MaxDiffK"],
        )
        for _, tract_row in tract_geodata.iterrows()
    }
    return subcounty_areas, tract_crosswalk


def save_snapshot(
    subcounty_areas: dict[str, SubcountyArea],
    tract_crosswalk: pd.DataFrame,
    layer_name: str,
    coordinate_reference_system,
) -> None:
    """Save both the current subcounty area map AND its crosswalk to original
    tracts, as two layers/tables in the same GeoPackage file."""
    snapshot_geodata = gpd.GeoDataFrame(
        [
            {
                "ClusterID": area.cluster_id,
                "ClusterPopulation": area.population,
                COUNTY_FIELD: area.county,
                "MaxDiffK": area.max_diff_k,
            }
            for area in subcounty_areas.values()
        ],
        geometry=[area.shape for area in subcounty_areas.values()],
        crs=coordinate_reference_system,
    )
    snapshot_geodata.to_file(OUTPUT_GEOPACKAGE_PATH, layer=layer_name, driver="GPKG", mode="w")

    database_connection = sqlite3.connect(OUTPUT_GEOPACKAGE_PATH)
    tract_crosswalk.to_sql(f"{layer_name}_Crosswalk", database_connection, if_exists="replace", index=False)
    database_connection.close()
    print(f"  Saved {layer_name} and its crosswalk table.")


def build_merge_diagnostic_record(
    county_mode: str,
    current_area_count: int,
    seed_id: str,
    seed_population: int,
    partner_id: str,
    chosen_candidate: MergeCandidate,
    ranked_all_pairs: list[MergeCandidate],
) -> dict:
    """One row of the merge-diagnostics CSV: the seed and its chosen partner,
    their population and clustering difference, the combined similarity score,
    and - critically - where that chosen match ranks (by percentile) against
    every pair available anywhere on the map that same round (including the
    chosen pair itself), so a "poor match" can be judged against what else was
    actually available that round rather than a fixed threshold."""
    all_population_diffs = np.array([p.population_difference for p in ranked_all_pairs])
    all_clustering_diffs = np.array([p.clustering_difference for p in ranked_all_pairs])
    all_scores = np.array([p.similarity_score for p in ranked_all_pairs])
    return {
        "county_mode": county_mode,
        "area_count_before_merge": current_area_count,
        "seed_id": seed_id,
        "seed_population": seed_population,
        "partner_id": partner_id,
        "population_difference": chosen_candidate.population_difference,
        "clustering_difference": chosen_candidate.clustering_difference,
        "similarity_score": chosen_candidate.similarity_score,
        "n_mapwide_candidates": len(ranked_all_pairs),
        "population_difference_percentile": float((all_population_diffs < chosen_candidate.population_difference).mean() * 100),
        "clustering_difference_percentile": float((all_clustering_diffs < chosen_candidate.clustering_difference).mean() * 100),
        "similarity_score_percentile": float((all_scores < chosen_candidate.similarity_score).mean() * 100),
    }


def combine_county_values(county_a: str, county_b: str) -> str:
    """Combine two areas' county values into one. Merging two same-county areas
    (the only kind COUNTY_CONSTRAINED_MODE ever merges) returns that single
    county unchanged; merging areas from different counties (possible only in
    COUNTY_UNCONSTRAINED_MODE) returns a sorted, deduplicated, comma-separated
    list of every county involved, rather than silently picking one and
    discarding the other."""
    combined_counties = set(county_a.split(",")) | set(county_b.split(","))
    return ",".join(sorted(combined_counties))


def apply_merge(
    subcounty_areas: dict[str, SubcountyArea],
    tract_crosswalk: pd.DataFrame,
    chosen_candidate: MergeCandidate,
    merged_shape: shapely.Geometry,
    merged_max_diff_k: float,
) -> None:
    """Fold chosen_candidate's two areas into one: point every original tract
    that belonged to either of them at the crosswalk to the merged cluster ID,
    remove the two input areas from subcounty_areas, and add the merged area
    in their place."""
    merged_cluster_id = min(chosen_candidate.area_id_a, chosen_candidate.area_id_b)
    area_a = subcounty_areas[chosen_candidate.area_id_a]
    area_b = subcounty_areas[chosen_candidate.area_id_b]
    merged_population = area_a.population + area_b.population
    merged_county = combine_county_values(area_a.county, area_b.county)

    tract_crosswalk.loc[
        tract_crosswalk["ClusterID"].isin([chosen_candidate.area_id_a, chosen_candidate.area_id_b]),
        "ClusterID",
    ] = merged_cluster_id

    del subcounty_areas[chosen_candidate.area_id_a]
    if chosen_candidate.area_id_b in subcounty_areas:
        del subcounty_areas[chosen_candidate.area_id_b]
    subcounty_areas[merged_cluster_id] = SubcountyArea(
        cluster_id=merged_cluster_id,
        shape=merged_shape,
        population=merged_population,
        county=merged_county,
        max_diff_k=merged_max_diff_k,
    )


def run_hierarchical_merge(
    tract_geodata: gpd.GeoDataFrame,
    county_mode: str,
    residential_sites: ResidentialSites,
) -> list[dict]:
    """Repeatedly merge the current lowest-population area into its best-matching
    touching neighbor until each target aggregation count is reached, maintaining
    a crosswalk to original tracts and a per-merge diagnostic record of how that
    match compared to everything else available on the map that round. Returns
    that round-by-round list of diagnostic records.

    If a round finds no eligible touching pairs before reaching the smallest
    target count - most likely in Counties mode, since nothing guarantees the
    remaining areas stay fully connected once merges are restricted to stay
    within a county, though NoCounties isn't structurally immune either (e.g.
    a genuinely disconnected landmass) - the current state is snapshotted
    anyway as the closest count actually reachable, rather than silently
    stopping without one."""
    coordinate_reference_system = tract_geodata.crs
    subcounty_areas, tract_crosswalk = build_initial_areas_and_crosswalk(tract_geodata)
    merge_diagnostics: list[dict] = []

    while True:
        current_area_count = len(subcounty_areas)
        is_target_count = current_area_count in AGGREGATION_TARGET_COUNTS

        if is_target_count:
            save_snapshot(
                subcounty_areas, tract_crosswalk, f"Subcounty{county_mode}{current_area_count}",
                coordinate_reference_system,
            )

        if current_area_count <= min(AGGREGATION_TARGET_COUNTS):
            break

        candidate_pairs = find_touching_area_pairs(subcounty_areas, county_mode)
        if not candidate_pairs:
            if not is_target_count:
                print(
                    f"  WARNING: no touching pairs remain ({county_mode}) at count "
                    f"{current_area_count}; saving this state as the closest count "
                    "reachable toward whichever smaller target(s) it never reached."
                )
                save_snapshot(
                    subcounty_areas, tract_crosswalk, f"Subcounty{county_mode}{current_area_count}",
                    coordinate_reference_system,
                )
            break

        ranked_all_pairs = rank_candidates_by_similarity(candidate_pairs)
        areas_with_candidate_pairs = {p.area_id_a for p in ranked_all_pairs} | {p.area_id_b for p in ranked_all_pairs}
        seed_id = select_lowest_population_seed(subcounty_areas, areas_with_candidate_pairs)

        seed_candidates = [p for p in ranked_all_pairs if p.area_id_a == seed_id or p.area_id_b == seed_id]
        chosen_candidate = seed_candidates[0]
        merged_shape = shapely.union(
            subcounty_areas[chosen_candidate.area_id_a].shape,
            subcounty_areas[chosen_candidate.area_id_b].shape,
        )

        # seed_id is always one side of chosen_candidate, so its population is already
        # known without re-checking which side; partner_id is whichever side isn't it.
        seed_population = subcounty_areas[seed_id].population
        partner_id = (
            chosen_candidate.area_id_b if chosen_candidate.area_id_a == seed_id
            else chosen_candidate.area_id_a
        )

        # Diagnostics: where does the chosen (seed-constrained) match rank against EVERY
        # pair available anywhere on the map this round - i.e. what a map-wide best-match
        # rule would have been choosing from instead.
        merge_diagnostics.append(
            build_merge_diagnostic_record(
                county_mode, current_area_count, seed_id, seed_population, partner_id,
                chosen_candidate, ranked_all_pairs,
            )
        )

        print(
            f"mode:{county_mode} | count:{current_area_count} | "
            f"seed:{seed_id} (pop {seed_population}) "
            f"merging {chosen_candidate.area_id_a} + {chosen_candidate.area_id_b}"
        )

        merged_max_diff_k = compute_max_diff_k(residential_sites, merged_shape)
        apply_merge(subcounty_areas, tract_crosswalk, chosen_candidate, merged_shape, merged_max_diff_k)

    return merge_diagnostics


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main() -> None:
    os.makedirs(os.path.dirname(OUTPUT_GEOPACKAGE_PATH), exist_ok=True)

    residential_sites, tract_geodata = prepare_datasets()

    print("Computing baseline Ripley's K for each tract...")
    baseline_max_diff_k_values = []
    for tract_polygon in tract_geodata.geometry:
        baseline_max_diff_k_values.append(compute_max_diff_k(residential_sites, tract_polygon))
    tract_geodata["MaxDiffK"] = baseline_max_diff_k_values

    print("Saving base tracts (before any merging) as their own snapshot...")
    base_areas, base_crosswalk = build_initial_areas_and_crosswalk(tract_geodata)
    save_snapshot(base_areas, base_crosswalk, "BaseTracts", tract_geodata.crs)

    merge_diagnostics: list[dict] = []
    for county_mode in COUNTY_MODES:
        merge_diagnostics.extend(run_hierarchical_merge(tract_geodata, county_mode, residential_sites))

    pd.DataFrame(merge_diagnostics).to_csv(MERGE_DIAGNOSTICS_CSV_PATH, index=False)
    print(f"Saved {len(merge_diagnostics)} merge diagnostic records to {MERGE_DIAGNOSTICS_CSV_PATH}")
    print(f"Done. Output written to {OUTPUT_GEOPACKAGE_PATH}")


if __name__ == "__main__":
    main()
