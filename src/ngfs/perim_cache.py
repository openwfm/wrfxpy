#caches perims vfrom arcgis NIFC portal
from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from datetime import datetime, timezone
import pandas as pd
from ingest.downloader import download_url
import geopandas as gpd
import simplekml
import hashlib
from pathlib import Path
import re
import requests


def download_ytd():
    file = None
    for sa in sys.argv:
        if '.geojson' in sa:
            file = sa
    if not file == None:
        gdf = load_geojson_with_processed_time(file)
        print(f'Using saved geojson file {file}')
        return gdf
    #downloads the year-to-date perimeter file
    url = 'https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/WFIGS_Interagency_Perimeters_YearToDate/FeatureServer/replicafilescache/WFIGS_Interagency_Perimeters_YearToDate_-237891591426750996.geojson'
          #https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/WFIGS_Interagency_Perimeters_YearToDate/FeatureServer/replicafilescache/WFIGS_Interagency_Perimeters_YearToDate_-237891591426750996.geojson
    now = pd.Timestamp.now('UTC')
    date_str = now.isoformat()[:10].replace(':','_') #like this '2026-03-10'
    save_str = f'ngfs/perims/perims_ytd_{date_str}.geojson'
    #download_url(url,save_str)
    changed = download_if_changed(url=url,local_path=save_str)
    if changed and os.path.exists(save_str):
        gdf = gpd.read_file(save_str)
        gdf = parse_dates(gdf)
        return gdf
    else:
        print('File not obtained')
        return gpd.GeoDataFrame()
    
    
def download_if_changed(url,local_path, timeout=30):
    """
    Download a GeoJSON file only if it differs from the local copy.
    Returns True if file was updated, False otherwise.
    """
    local_path = Path(local_path)
    # Fetch remote content
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    remote_content = response.content
    # Compute remote hash
    remote_hash = hashlib.sha256(remote_content).hexdigest()
    # Compare with local file (if it exists)
    if local_path.exists():
        local_hash = file_hash(local_path)
        if local_hash == remote_hash:
            print("No change detected. Skipping download.")
            return False
    # Write new file
    local_path.parent.mkdir(parents=True, exist_ok=True)
    with open(local_path, "wb") as f:
        f.write(remote_content)
    print("File updated.")
    return True
    
def find_incident_files(repo_dir = '/data/jhaley/wrfxpy/ngfs/perims', irwin_id = '123123123123'):
    """
    Return a list of all GeoJSON files matching the IrwinID.
    """
    repo_dir = Path(repo_dir)
    matches = list(repo_dir.glob(f"*_{irwin_id}.geojson"))
    return matches
    
def make_filename(incident_name, irwin_id):
    """
    Create human-readable, filesystem-safe filename.
    """
    # Remove unsafe characters
    safe_name = re.sub(r'[^A-Za-z0-9_\-]+', '_', incident_name.strip())
    return f"{safe_name}_{irwin_id}.geojson"

def parse_dates(gdf):
    '''
    finds date strings in GeoDataframe gdf and vonverts to time object
    '''
    for k in gdf.keys():
        if 'time' in k.lower() or 'date' in k.lower() or k == 'processed_utc':
            gdf[k] = pd.to_datetime(gdf[k], utc=True, errors="coerce")
    return gdf

def auto_parse_datetime_columns(gdf, utc=True, threshold=0.8):
    """
    Automatically convert timestamp-like columns in a GeoDataFrame
    to pandas datetime dtype.
    Parameters
    ----------
    gdf : GeoDataFrame
    utc : bool
        Convert to UTC timezone.
    threshold : float
        Fraction of values that must successfully parse to treat
        the column as datetime.
    """
    gdf = gdf.copy()
    for col in gdf.columns:
        # Skip geometry
        if col == gdf.geometry.name:
            continue
        # Only attempt parsing object/string columns
        if gdf[col].dtype != "object":
            continue
        parsed = pd.to_datetime(gdf[col], utc=utc, errors="coerce")
        success_fraction = parsed.notna().mean()
        if success_fraction >= threshold:
            gdf[col] = parsed
    return gdf
        

def write_kml_polygon_with_attributes(gdf, out_path):

    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")

    kml = simplekml.Kml()

    for _, row in gdf.iterrows():
        geom = row.geometry

        if geom.geom_type == "Polygon":
            pol = kml.newpolygon()
            pol.outerboundaryis = list(geom.exterior.coords)

        elif geom.geom_type == "Point":
            pol = kml.newpoint()
            pol.coords = [(geom.x, geom.y)]

        else:
            continue  # extend as needed

        # Add all non-geometry fields
        for col in gdf.columns:
            if col != "geometry":
                value = row[col]
                if value is not None:
                    pol.extendeddata.newdata(name=col, value=str(value))

    kml.save(str(out_path))

def write_kml_with_attributes(gdf, out_path):
    print(f'Writing kml file {out_path}')
    # Ensure KML-compatible CRS
    #if gdf.crs != "EPSG:4326":
    #    gdf = gdf.to_crs("EPSG:4326")

    kml = simplekml.Kml()

    for _, row in gdf.iterrows():
        geom = row.geometry

        # ---- POLYGON ----
        if geom.geom_type == "Polygon":
            polys = [geom]

        # ---- MULTIPOLYGON ----
        elif geom.geom_type == "MultiPolygon":
            polys = list(geom.geoms)
        else:
            continue  # Extend if you need LineString/Point
        for poly in polys:
            pol = kml.newpolygon()
            pol.outerboundaryis = list(poly.exterior.coords)
            # Add attributes
            for col in gdf.columns:
                if col != "geometry":
                    value = row[col]
                    if value is not None:
                        pol.extendeddata.newdata(
                            name=col,
                            value=str(value)
                        )
    kml.save(str(out_path))

def file_hash(path):
    """Compute SHA256 hash of a local file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def geometry_hash(geom):
    """
    Return a stable hash for a shapely geometry.
    """
    wkb = geom.wkb  # binary representation
    return hashlib.sha256(wkb).hexdigest()
    
def update_incident_store_geojson(gdf,
                                  repo_dir = '/data/jhaley/wrfxpy/ngfs/perims',
                                  irwin_col="poly_IRWINID",
                                  name_col = 'attr_IncidentName',
                                  id_col="OBJECTID"):
    repo_dir = Path(repo_dir)
    repo_dir.mkdir(parents=True, exist_ok=True)
    # Add reliable UTC timestamp
    gdf = gdf.copy()
    if 'processed_utc' not in gdf.keys():
        gdf["processed_utc"] = datetime.now(timezone.utc)
    print(f'Geodataframe has {len(gdf)} lines to process')
    for irwin_id, group in gdf.groupby(irwin_col):
        #try to get the ID from attr_IrwinID
        if pd.isna(irwin_id):
            irwin_id = group['attr_IrwinID']
            if pd.isna(irwin_id):
                irwin_id = "UNASSIGNED"
        out_files = find_incident_files(repo_dir=repo_dir,irwin_id=irwin_id)
        if len(out_files) > 1:
            print(f"Note: {len(out_files)} files found for IrwinID {irwin_id}")
        if len(out_files):
            for out_file in out_files:
                existing = gpd.read_file(out_file)
                existing = parse_dates(gdf)
                # Remove already stored OBJECTIDs
                new_rows = group[~group[id_col].isin(existing[id_col])]
                if len(new_rows) > 0:
                    print(f'Updating {out_file}')
                    combined = pd.concat([existing, new_rows], ignore_index=True)
                    print(combined.processed_utc.unique())
                    combined = combined.sort_values("processed_utc")
                    combined.to_file(out_file, driver="GeoJSON")
                    write_kml_with_attributes(combined,out_file.replace('.geojson','.kml'))
        else:
            out_name = make_filename(incident_name=group[name_col].iloc[0],irwin_id=irwin_id)
            out_file = f'{repo_dir}/{out_name}'
            print(f'Saving {out_file}')
            group.to_file(out_file, driver="GeoJSON")
            kml_file = out_file.replace('.geojson','.kml')
            try:
                write_kml_with_attributes(group,kml_file)
            except:
                print(f'Error writing kml file {kml_file}')
                print(group)
    
def save_if_geometry_changed(gdf, repo_dir,
                             id_col="OBJECTID"):
    repo_dir = Path(repo_dir)
    repo_dir.mkdir(parents=True, exist_ok=True)

    saved_paths = []

    for _, row in gdf.iterrows():
        obj_id = str(row[id_col])
        geom = row.geometry

        geom_hash = geometry_hash(geom)

        obj_dir = repo_dir / obj_id
        obj_dir.mkdir(exist_ok=True)

        out_file = obj_dir / f"{geom_hash}.gpkg"

        if not out_file.exists():
            single = gdf.loc[[row.name]].copy()
            single["processed_utc"] = row.get("processed_utc")
            single.to_file(out_file, driver="GPKG")

        saved_paths.append(out_file)

    return saved_paths
    

def load_geojson_with_processed_time(geojson_path):
    """
    Load a GeoJSON file and add a 'processed_utc' column using the
    file's last modification time (UTC).
    """
    geojson_path = Path(geojson_path)
    # Read file modification time
    mtime = geojson_path.stat().st_mtime
    processed_utc = datetime.fromtimestamp(mtime, tz=timezone.utc)
    # Load GeoJSON
    gdf = gpd.read_file(geojson_path)
    gdf = parse_dates(gdf)
    # Add column
    if 'processed_utc' not in gdf.keys():
        gdf["processed_utc"] = processed_utc
        print(f'Adding modified time {processed_utc}')
    return gdf


if __name__ == "__main__":
    print(f"Starting perimeter caching {pd.Timestamp.now('UTC')}")
    gdf = download_ytd()
    if len(gdf) > 0:
        print('Processing new data')
        update_incident_store_geojson(gdf)
    else:
        print('No new data available')
    print(f"Finished perimeter caching {pd.Timestamp.now('UTC')}")