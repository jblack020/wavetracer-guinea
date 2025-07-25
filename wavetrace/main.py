"""
CONVENTIONS:
    - All longitudes and latitudes below are referenced to the WGS84 ellipsoid, unless stated otherwise
"""
from pathlib import Path
import re
import csv
import textwrap
import shutil
import subprocess
import base64
from math import sin, cos, atan, atan2, sqrt, pi, radians, degrees
import tempfile
import urllib
import time
from shapely.geometry import Point
import requests
from tqdm import tqdm
import wavetrace.constants as cs
import wavetrace.utilities as ut
import pandas as pd
import geopandas as gpd
import numpy as np
import subprocess
import tempfile
import shutil
import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio
from rasterio import features
from shapely.geometry import shape
from pathlib import Path


# In the calls to the subprocess function below,
# sometimes instead of using absolute paths,
# i use relative paths in combination with the ``cwd`` option.
# Mind the difference!

def process_transmitters(in_path, out_path,
                         earth_dielectric_constant=cs.EARTH_DIELECTRIC_CONSTANT,
                         earth_conductivity=cs.EARTH_CONDUCTIVITY,
                         radio_climate=cs.RADIO_CLIMATE,
                         fraction_of_time=cs.FRACTION_OF_TIME,
                         fraction_of_situations=cs.FRACTION_OF_SITUATIONS):
    """
    Read the CSV transmitter data at ``in_path``, and for each transmitter, 
    create the following SPLAT! data for it and save it to the directory
    ``out_path``:

    - location data as a QTH file
    - irregular topography model parameter as an LRP file
    - azimuth data as an AZ file
    - elevation data as an EL file

    INPUT:
        - ``in_path``: string or Path object; location of a CSV file of transmitter data
        - ``out_path``: string or Path object; directory to which to write outputs
        - ``earth_dielectric_constant``: float; Earth dielectric constant; SPLAT! parameter used to make an LRP file
        - ``earth_conductivity``: float; Earth conductivity; SPLAT! parameter used to make an LRP file
        - ``radio_climate``: integer; SPLAT! parameter used to make an LRP file
        - ``fraction_of_time``: float in [0, 1]; SPLAT! parameter used to make an LRP file
        - ``fraction_of_situations``: float in [0, 1]; SPLAT! parameter used to make an LRP file

    OUTPUT:
        None.

    NOTES:
        The CSV file of transmitter data must include at least the columns,
        otherwise a ``ValueError`` will be raised.

        - ``'network_name'``: name of transmitter network
        - ``'site_name'``: name of transmitter site
        - ``'longitude'``: decimal longitude of transmitter  
        - ``'latitude``: decimal latitude of transmitter
        - ``'antenna_height'``: height of transmitter antenna in meters above sea level
        - ``'polarization'``: 0 for horizontal or 1 for vertical
        - ``'frequency'``: frequency of transmitter in megaherz
        - ``'power_eirp'``: effective radiated power of transmitter in watts
    """
    # Read transmitter data
    ts = read_transmitters(in_path)

    # Write SPLAT files
    out_path = Path(out_path)
    if not out_path.exists():
        out_path.mkdir(parents=True)

    for t in ts:
        for f, kwargs, ext in [
            (build_splat_qth, {}, '.qth'),
            (build_splat_lrp, {
                'earth_dielectric_constant': earth_dielectric_constant,
                'earth_conductivity': earth_conductivity,
                'radio_climate': radio_climate,
                'fraction_of_time': fraction_of_time
            }, '.lrp'),
            (build_splat_az, {}, '.az'),
            (build_splat_el, {}, '.el'),
        ]:
            s = f(t, **kwargs)
            path = out_path/(t['name'] + ext)
            with path.open('w') as tgt:
                tgt.write(s)


def read_transmitters(path):
    """
    Return a list of dictionaries, one for each transmitter in the transmitters
    CSV file.

    INPUT:
        - ``path``: string or Path object; location of a CSV file of transmitter data

    OUTPUT:
        List of dictionaries.
        The keys for each transmitter come from the header row of the CSV file.
        If ``REQUIRED_TRANSMITTER_FIELDS`` is not a subset of these keys, then
        raise a ``ValueError``.
        Additionally, a 'name' field is added to each transmitter dictionary for later use and is the result of :func:`build_transmitter_name`.

    NOTES:
        For the format of the transmitters CSV file, see the notes section of :func:`create_splat_transmitters`.
    """
    path = Path(path)
    transmitters = []
    with path.open() as src:
        reader = csv.DictReader(src)
        for row in reader:
            transmitters.append(row)
    transmitters = check_and_format_transmitters(transmitters)
    return transmitters


def check_and_format_transmitters(transmitters):
    """
    Check and format the given list of transmitter dictionaries.

    INPUT:
        - ``transmitters``: list; same format as output of :func:`read_transmitters`

    OUTPUT:
        The given list of transmitters dictionaries altered as follows.
        For each dictionary, 

        - create a ``name`` field from the ``network_name`` and ``site_name`` fields
        - convert the numerical fields to floats

    NOTES:
        Raises a ``ValueError`` if the list of transmitters is empty, or if the ``REQUIRED_TRANSMITTER_FIELDS`` are not present in each transmitter dictionary, or if the any of the field data is improperly formatted.
    """
    if not transmitters:
        raise ValueError('Transmitters must be a nonempty list')

    # Check that required fields are present
    keys = transmitters[0].keys()
    if not set(cs.REQUIRED_TRANSMITTER_FIELDS) <= set(keys):
        raise ValueError('Transmitters header must contain '
                         'at least the fields {!s}'.format(cs.REQUIRED_TRANSMITTER_FIELDS))

    # Format required fields and raise error if run into problems
    new_transmitters = []
    for i, t in enumerate(transmitters):
        try:
            new_name = build_transmitter_name(t['network_name'],
                                              t['site_name'])
            i = 1
            while new_name in [t['name'] for t in new_transmitters]:
                new_name = new_name + '_' + str(i)
                i += 1
            t['name'] = new_name

            for key in ['latitude', 'longitude', 'antenna_height',
                        'polarization', 'frequency', 'power_eirp']:
                t[key] = float(t[key])
        except:
            raise ValueError('Data on line {!s} of transmitters file is '
                             'improperly formatted'.format(i + 1))
        new_transmitters.append(t)

    return new_transmitters


def build_transmitter_name(network_name, site_name):
    """
    Return a string that is the network name with spaces removed followed 
    by an underscore followed by the site name with spaces removed.

    EXAMPLES:

    >>> build_transmitter_name('Slap hAppy', 'Go go ')
    'SlaphAppy_Gogo'
    """
    return network_name.replace(' ', '') + '_' +\
        site_name.replace(' ', '')


def build_splat_qth(transmitter):
    """
    Return the text content of a SPLAT! site location file 
    (QTH file) corresponding to the given transmitter.

    INPUT:
        - ``transmitter``: dictionary of the same form as any one of the elements in the list output by :func:`read_transmitters`

    OUTPUT:
        String.
    """
    t = transmitter
    # Convert to degrees east in range (-360, 0] for SPLAT!
    lon = -t['longitude']
    return "{!s}\n{!s}\n{!s}\n{!s}m".format(
        t['name'],
        t['latitude'],
        lon,
        t['antenna_height'])


def build_splat_lrp(transmitter,
                    earth_dielectric_constant=cs.EARTH_DIELECTRIC_CONSTANT,
                    earth_conductivity=cs.EARTH_CONDUCTIVITY,
                    radio_climate=cs.RADIO_CLIMATE,
                    fraction_of_time=cs.FRACTION_OF_TIME,
                    fraction_of_situations=cs.FRACTION_OF_SITUATIONS):
    """
    Return the text (string) content of a SPLAT! irregular topography model parameter file (LRP file) corresponding to the given transmitter.

    INPUT:
        - ``transmitter``: dictionary of the same form as any one of the elements in the list output by :func:`read_transmitters`
        - ``earth_dielectric_constant``: float
        - ``earth_conductivity``: float
        - ``radio_climate``: integer
        - ``fraction_of_time``: float in [0, 1]
        - ``fraction_of_situations``: float in [0, 1]

    OUTPUT:
        String
    """
    t = transmitter
    s = """\
    {!s} ; Earth Dielectric Constant (Relative permittivity)
    {!s} ; Earth Conductivity (Siemens per meter)
    301.000 ; Atmospheric Bending Constant (N-units)
    {!s} ; Frequency in MHz (20 MHz to 20 GHz)
    {!s} ; Radio Climate
    {!s} ; Polarization (0 = Horizontal, 1 = Vertical)
    {!s} ; Fraction of situations
    {!s} ; Fraction of time 
    {!s} ; ERP in watts""".format(
        earth_dielectric_constant,
        earth_conductivity,
        t['frequency'],
        radio_climate,
        t['polarization'],
        fraction_of_situations,
        fraction_of_time,
        t['power_eirp'])
    return textwrap.dedent(s)


def build_splat_az(transmitter):
    """
    Return the text (string) content of a SPLAT! azimuth file (AZ file)
    corresponding to the given transmitter.

    If both ``bearing`` (boresight direction) and ``horizontal_beamwidth`` are
    present, we create a simple three-level pattern:

        • 0 dB (linear = 1.0) inside the main beam  
        • 14 dB (linear ≈ 0.2) in the forward side-lobes  
        • 26 dB (linear ≈ 0.05) in the rear lobe  

    Otherwise we fall back to the SPLAT! placeholder ``0  0`` (omnidirectional).
    """
    try:
        bearing = float(transmitter["bearing"]) % 360          # degrees
        hb = float(transmitter["horizontal_beamwidth"])    # degrees

        MAIN_GAIN = 1.0     # 0 dB
        SIDE_GAIN = 0.2     # –14 dB
        BACK_GAIN = 0.05    # –26 dB
        half_bw = hb / 2

        # first line = boresight
        lines = [f"{bearing}"]
        for deg in range(360):
            # angle from boresight
            rel = (deg - bearing) % 360
            if rel <= half_bw or rel >= 360 - half_bw:          # main lobe
                g = MAIN_GAIN
            elif rel < 180:                                     # forward side-lobes
                g = SIDE_GAIN
            else:                                               # rear lobe
                g = BACK_GAIN
            lines.append(f"{deg}  {g:.3f}")
        return "\n".join(lines)
    except (KeyError, ValueError, TypeError):
        # Missing or bad pattern info → omnidirectional placeholder
        return "0  0"


def build_splat_el(transmitter):
    """
    Return the text (string) content of a SPLAT! elevation file (EL file) corresponding to the given transmitter.

    INPUT:
        - ``transmitter``: dictionary of the same form as any one of the elements in the list output by :func:`read_transmitters`

    OUTPUT:
        String

    NOTES:
        A transmitter with no ``'bearing'`` or ``'antenna_downtilt'`` or ``'vertical_beamwidth'`` data will produce the string ``'0  0'``.
    """
    t = transmitter
    try:
        bearing = float(t['bearing'])
        ad = float(t['antenna_downtilt'])
        vb = float(t['vertical_beamwidth'])
        s = '{!s}  {!s}\n'.format(ad, bearing)
        counter = 0
        for x in range(-10, 91):
            if counter < vb:
                s += '{!s}  0.9\n'.format(x)
            else:
                s += '{!s}  0.1\n'.format(x)
            counter += 1
    except:
        s = '0  0\n'

    return s[:-1]  # Drop the final newline


def get_lonlats(transmitters):
    """
    Return a list of longitude-latitude pairs (float pairs) representing the locations of the given transmitters.
    If ``transmitters`` is empty, then return the empty list. 

    INPUT:
        - ``transmitters``: a list of transmitters of the form output by :func:`read_transmitters`

    OUTPUT:
        String
    """
    return [(t['longitude'], t['latitude']) for t in transmitters]


def get_covering_tiles_ids(transmitters, transmitter_buffer=0.5):
    """
    Given a list of transmitters (of the form output by :func:`read_transmitters`), get their locations, buffer them by ``transmitter_buffer`` decimal degrees, and return an ordered list of the unique SRTM tile IDs in ``tile_ids`` whose corresponding tiles intersect the buffers.
    As long as ``tile_ids`` and ``transmitter_buffer`` are big enough, the result will be a list of tile IDs to use when computing coverage for the given transmitters.
    The defaults are appropriate for transmitters in Guinea.

    NOTES:
        - Regarding the transmitter buffer, one degree of latitude represents about 111 km on the ground and one degree of longitude at -45 degrees latitude represents about 78 km on the ground; see https://en.wikipedia.org/wiki/Decimal_degrees
    """
    blobs = [Point(p).buffer(transmitter_buffer)
             for p in get_lonlats(transmitters)]
    return ut.compute_intersecting_tiles(blobs)


def download_topography(tile_ids, path, high_definition=False):
    """
    Download from the public Gitlab repository https://gitlab.com/jblack020-group/srtm_data_guinea the SRTM3 topography data corresponding to the given SRTM tile IDs and save the files to the directory ``path``, creating the directory if it does not exist.

    INPUT:
        - ``tile_ids``: list of strings; SRTM tile IDs
        - ``path``: string or Path object specifying a directory
        - ``high_definition``: boolean; if ``True`` then download SRTM1 tiles; otherwise download SRTM3 tiles

    OUTPUT:
        None

    NOTES:
        Only works for SRTM tiles covering Guinea and raises a ``ValueError`` if the set of tile IDs is not a subset of :data:`SRTM_GIN_TILE_IDS`
    """
    if not set(tile_ids) <= set(cs.SRTM_GIN_TILE_IDS):
        raise ValueError("Tile IDs must be a subset of {!s}".format(
            ' '.join(cs.SRTM_GIN_TILE_IDS)))

    # Set download parameters
    project_id = '70958368'
    url = 'https://gitlab.com/api/v4/projects/{!s}/repository/files/'.\
        format(project_id)
    if high_definition:
        file_names = ['srtm1/{!s}.SRTMGL1.hgt.zip'.format(t) for t in tile_ids]
    else:
        file_names = ['srtm3/{!s}.SRTMGL3.hgt.zip'.format(t) for t in tile_ids]

    # Create output directory
    path = Path(path)
    if not path.exists():
        path.mkdir(parents=True)

    # Download
    params = {'ref': 'main', }
    for file_name in tqdm(file_names, total=len(file_names), desc="Downloading topography data"):
        time.sleep(20)
        file_url = '{url}{file_id}'.format(
            url=url,
            file_id=urllib.parse.quote_plus(file_name))
        r = requests.get(file_url, params=params, stream=True)

        if r.status_code != requests.codes.ok:
            raise ValueError('Downloading file {!s} failed with status '
                             ' code {!s}'.format(file_name, r.status_code))

        p = path/file_name.split('/')[-1]
        with p.open('wb') as tgt:
            content = base64.b64decode(r.json()['content'])
            tgt.write(content)


def process_topography(in_path, out_path, high_definition=False):
    """
    Convert each SRTM HGT topography file in the directory ``in_path`` to a SPLAT! Data File (SDF) file in the directory ``out_path``,     creating the directory if it does not exist.
    If ``high_definition``, then assume the input data is high definition.

    INPUT:
        - ``in_path``: string or Path object specifying a directory
        - ``out_path``: string or Path object specifying a directory
        - ``high_definition``: boolean

    OUTPUT:
        None.

    NOTES:
        - Calls SPLAT!'s ``srtm2sdf`` or ``srtm2sdf-hd`` 
          (if ``high_definition``) command to do the work
        - Raises a ``subprocess.CalledProcessError`` if SPLAT! fails to 
          convert a file
        - Each SRTM1 or SRTM3 file must have a name of the form <SRTM tile ID>[.something].hgt.zip or <SRTM tile ID>[.something].hgt, e.g. S36E173.SRTMGL3.hgt.zip 
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    if not out_path.exists():
        out_path.mkdir(parents=True)

    splat = 'srtm2sdf'
    if high_definition:
        splat += '-hd'

    sdf_pattern = re.compile(r"[\d\w\-\:]+\.sdf")

    for f in tqdm(in_path.iterdir(), total=len(list(in_path.iterdir())), desc="Converting topography data to SDF"):
        if not (f.name.endswith('.hgt') or f.name.endswith('.hgt.zip')):
            continue

        # Unzip if necessary
        is_zip = False
        if f.name.endswith('.zip'):
            is_zip = True
            shutil.unpack_archive(str(f), str(f.parent))
            tile_id = f.name.split('.')[0]
            f = f.parent/'{!s}.hgt'.format(tile_id)

        # Convert to SDF
        cp = subprocess.run([splat, f.name], cwd=str(f.parent),
                            stdout=subprocess.PIPE, universal_newlines=True, check=True)

        # Get name of output file, which SPLAT! created and which differs
        # from the original name, and move the output to the out path
        m = sdf_pattern.search(cp.stdout)
        name = m.group(0)
        src = in_path/name
        tgt = out_path/name
        shutil.move(str(src), str(tgt))

        # Clean up
        if is_zip:
            f.unlink()

# @ut.time_it


def compute_coverage_0(in_path, out_path, transmitters,
                       high_definition=False):
    """
    Create a SPLAT! coverage report for every transmitter with data located at ``in_path``, or if ``transmitters`` is given, then every transmitter 
    in that list with data data located at ``in_path``.
    Write each report to the directory ``out_path``, creating the directory   if necessary.
    A report comprises the files

    - ``'<transmitter name>-site_report.txt'``
    - ``'<transmitter name>.kml'``: KML file containing transmitter feature and ``'<transmitter name>.ppm'``
    - ``'<transmitter name>.ppm'``: PPM file depicting a contour plot of the transmitter signal strength
    - ``'<transmitter name>-ck.ppm'``: PPM file depicting a legend for the signal strengths in ``'<transmitter name>.ppm'``

    INPUT:
        - ``in_path``: string or Path object specifying a directory; all the SPLAT! transmitter and elevation data should lie here
        - ``out_path``: string or Path object specifying a directory
        - ``transmitters``: list of transmitter dictionaries (in the form output by :func:`read_transmitters`) to grab the frequencies of to convert dBμV/m to dBm (SPLAT! uses dBm)
        - ``receiver_sensitivity``: float; desired path loss threshold beyond which signal strength contours will not be plotted (measured in dBμV/m)
        - ``high_definition``: boolean

    OUTPUT:
        None. 

    NOTES:
        - Calls SPLAT!'s ``splat`` or ``splat-hd`` (if ``high_definition``) to do the work
        - Raises a ``subprocess.CalledProcessError`` if SPLAT! fails
        - This is a time-intensive function. 
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    if not out_path.exists():
        out_path.mkdir(parents=True)

    tx_index = {t['name']: t for t in transmitters or []}

    # Get transmitter names
    if transmitters is not None:
        transmitter_names = [t['name'] for t in transmitters]
    else:
        transmitter_names = [p.stem for p in in_path.iterdir()
                             if p.name.endswith('.qth')]

    print(f"Transmitter names: {transmitter_names}")

    # Splatify
    splat = 'splat'
    if high_definition:
        splat += '-hd'

    for t in tqdm(transmitter_names, total=len(transmitter_names), desc="Computing coverage"):
        power_erp = tx_index[t]['power_eirp']  # ERP

        # build splat argument list
        args = [
            splat,
            '-t', f'{t}.qth',
            '-L', '8.0',  # Rx height (m or ft)
            '-metric', '-ngs',  # Use metres, hide greyscale topo
            '-kml',  # Google-Earth overlay (optional)
            '-ano', f'{t}.dat',  # alphanumeric output file
            '-o', f'{t}.ppm',  # Colour contour map
            '-erp', str(power_erp)  # ERP
            
        ]

        subprocess.run(
            args,
            cwd=str(in_path),
            stdout=subprocess.PIPE,
            universal_newlines=True,
            check=True
        )

    # Move outputs to out_path
    exts = ['.ppm', '-ck.ppm', '-site_report.txt', '.kml', '.dat']
    for t in transmitter_names:
        for ext in exts:
            src = in_path/(t + ext)
            tgt = out_path/(t + ext)
            try:
                shutil.move(str(src), str(tgt))
            except FileNotFoundError:
                print(f"File {src} not found")


def postprocess_coverage_0(path, keep_ppm, make_shp):
    """
    Using the PPM files in the directory ``path`` do the following:

    - Convert each PPM file into a PNG file, replacing white with transparency using ImageMagick
    - Change the PPM reference in each KML file to the corresponding PNG file
    - Convert the PNG coverage file (not the legend file) into GeoTIFF using GDAL
    - Optionally create ESRI Shapefile bundles (.dbf, .prj, .shp, and .shx files) from the GeoTIFF files

    INPUT:
        - ``path``: string or Path object; directory where coverage reports (outputs of :func:`compute_coverage`) lie
        - ``keep_ppm``: boolean; keep the original, large PPM files in the coverage reports if and only if this flag is ``True``
        - ``make_shp``: boolean; create ESRI Shapefiles from the GeoTIFF files if and only if this flag is ``True``

    OUTPUT:
        None.
    """
    path = Path(path)

    # First pass: create PNG from PPM
    for f in path.iterdir():
        if f.suffix == '.ppm':
            # Convert to PNG, turning white background into
            # transparent background
            png = f.stem + '.png'
            args = ['convert', '-transparent', '#FFFFFF', f.name, png]
            subprocess.run(args, cwd=str(path),
                           stdout=subprocess.PIPE, universal_newlines=True, check=True)

            # # Resize to width 1200 pixels
            # args = ['convert', '-geometry', '1200', png, png]
            # subprocess.run(args, cwd=str(path),
            #   stdout=subprocess.PIPE, universal_newlines=True, check=True)

            if not keep_ppm:
                # Delete PPM
                f.unlink()

    # Second pass: create KML and convert PNG to GeoTIFF
    for f in path.iterdir():
        if f.suffix == '.kml':
            # Replace PPM with PNG in KML
            with f.open() as src:
                kml = src.read()
            kml = kml.replace('.ppm', '.png')
            with f.open('w') as tgt:
                tgt.write(kml)

            # Convert main PNG to GeoTIFF using the lon-lat bounds from the KML
            bounds = get_bounds_from_kml(kml)
            ulx, uly, lrx, lry = str(bounds[0]), str(bounds[3]), \
                str(bounds[2]), str(bounds[1])
            epsg = 'EPSG:4326'  # WGS84
            png = f.stem + '.png'
            tif = f.stem + '.tif'
            args = ['gdal_translate', '-of', 'Gtiff', '-a_ullr',
                    ulx, uly, lrx, lry, '-a_srs', epsg, png, tif]
            subprocess.run(args, cwd=str(path),
                           stdout=subprocess.PIPE, universal_newlines=True, check=True)

    # Optional third pass: create vector files from GeoTIFFs
    if make_shp:
        for f in path.iterdir():
            if f.suffix == '.tif':
                tif = f.name
                shp = f.stem + '.shp'
                args = ['gdal_polygonize.py', tif, '-f', 'ESRI Shapefile', shp]
                subprocess.run(args, cwd=str(path),
                               stdout=subprocess.PIPE, universal_newlines=True, check=True)


def post_process_coverage_1(dat_path,
                            receiver_sensitivity,
                            out_dir=None,
                            crs="EPSG:4326",
                            pixel_deg=0.0025):
    """
    Read SPLAT *.dat* and keep val > threshold. Burn those cells into a boolean raster. Polygonise contiguous True cells. Write Shapefile + PNG.
    """
    dat_path = Path(dat_path)
    if out_dir is None:
        out_dir = dat_path.parent
    out_dir = Path(str(out_dir) + "_vector")
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = dat_path.stem

    # ── 1. read + filter ────────────────────────────────────────────────
    cols = ["lat_n", "lon_w", "az", "el", "val"]
    df = (pd.read_csv(dat_path, skiprows=2, names=cols)
          .apply(pd.to_numeric, errors="coerce")
          .dropna())
    df["lon"] = -df.lon_w
    df = df[(df.el > -90) & (df.val > receiver_sensitivity)]
    if df.empty:
        raise ValueError("Threshold filters out every point.")

    # ── 2. build target raster grid ─────────────────────────────────────
    xmin, xmax = df.lon.min(), df.lon.max()
    ymin, ymax = df.lat_n.min(), df.lat_n.max()
    width = int(np.ceil((xmax - xmin) / pixel_deg))
    height = int(np.ceil((ymax - ymin) / pixel_deg))
    transform = rasterio.transform.from_origin(
        xmin, ymax, pixel_deg, pixel_deg)

    # ── 3. burn points ----------------------------------------------------------------
    # Snap each point to its row / col index; duplicates are cheap
    rows = ((ymax - df.lat_n) / pixel_deg).astype(int)
    cols = ((df.lon - xmin) / pixel_deg).astype(int)
    raster = np.zeros((height, width), dtype="uint8")
    raster[rows, cols] = 1          # set covered cells to 1

    # ── 4. polygonise contiguous cells -----------------------------------------------
    shapes = features.shapes(
        raster, mask=raster.astype(bool), transform=transform)
    polys = [shape(geom) for geom, value in shapes if value == 1]

    gdf = gpd.GeoDataFrame(geometry=polys, crs=crs)
    shp = out_dir / f"{stem}.shp"
    gdf.to_file(shp)
    print(f"✓ wrote {shp} ({len(gdf)} polygons)")

    # optional PNG preview
    import imageio
    import matplotlib.pyplot as plt
    png = out_dir / f"{stem}.png"
    imageio.imwrite(png, (raster[::-1] * 255).astype("uint8"))
    print(f"✓ wrote {png}")


def get_bounds_from_kml(kml_string):
    """
    Given the text content of a SPLAT! KML coverage file, return a list of floats of the form ``[min_lon, min_lat, max_lon, max_lat]`` which describes the longitude-latitude bounding box of the coverage file.
    Raise an ``AttributeError`` if the KML does not contain a ``<LatLonBox>``  entry and hence is not a well-formed SPLAT! KML coverage file.
    """
    kml = kml_string
    west = re.search(r"<west>([0-9-][0-9\.]*)<\/west>", kml).group(1)
    south = re.search(r"<south>([0-9-][0-9\.]*)<\/south>", kml).group(1)
    east = re.search(r"<east>([0-9-][0-9\.]*)<\/east>", kml).group(1)
    north = re.search(r"<north>([0-9-][0-9\.]*)<\/north>", kml).group(1)
    result = [west, south, east, north]
    return list(map(float, result))


def compute_coverage(in_path, out_path, transmitters=None,
                     receiver_sensitivity=cs.RECEIVER_SENSITIVITY,
                     high_definition=False):
    """
    Produce coverage reports by running :func:`compute_coverage_0` and then run :func:`post_process_coverage_0`.
    """
    compute_coverage_0(in_path, out_path, transmitters,
                       high_definition)

    # Process each .dat file individually
    out_path = Path(out_path)
    for dat_file in out_path.glob("*.dat"):
        post_process_coverage_1(dat_file, receiver_sensitivity)


def partition(width, height, n=3):
    """
    Given the pixel width and pixel height of a rectangular image and an integer ``n``, partition the rectangle into ``n**2`` subrectangles, each of roughly the same sizes.
    Return a list of the subrectangle offsets and sizes for easy use with  GDAL's ``gdal_translate -srcwin`` option.
    Each list item has the form (x-offset, y-offset, x-size, y-size) and the items/subrectangles are ordered from left to right and then from top to bottom, e.g. for ``n=3`` the layout of the subrectangles looks like this::

        -------------
        | 0 | 1 | 2 |
        -------------
        | 3 | 4 | 5 |
        -------------
        | 6 | 7 | 8 |
        -------------

    The subrectangles in the right-most column and those in the bottom-most row will be slightly wider and taller, respectively, than the other subrectangles in case ``width`` or ``height`` are not divisible by ``n``, respectively.
    """
    q, r = divmod(width, n)
    xs = [(i*q, q) for i in range(n - 1)] + [((n - 1)*q, q + r)]
    q, r = divmod(height, n)
    ys = [(i*q, q) for i in range(n - 1)] + [((n - 1)*q, q + r)]
    return [(xoff, yoff, xsize, ysize) for yoff, ysize in ys
            for xoff, xsize in xs]


def compute_look_angles(lon, lat, height, satellite_lon):
    """
    Given the longitude, latitude, and height in meters of a point P on Earth and given the longitude of a geostationary satellite S, return the azimuth and elevation in degrees of S relative to P, respectively.

    INPUT:
        - ``lon``: float; longitude of P 
        - ``lat```: float; latitude of P
        - ``height``: float; distance in meters between P and the WGS84 ellipsoid; GPS height
        - ``satellite_lon``: float; longitude of S

    OUTPUT:
        - ``azimuth``: float; degrees in [0, 360)
        - ``elevation``: float; degrees in [-90, 90]; a negative value indicates that S lies below the local horizon of P

    NOTES:

    - Algorithm taken from `Determination of look angles to geostationary communication satellites <https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf>`_ by Tomas Soler David W. Eisemann
    - The input ``height`` is the sum of H and N, where H is the SRTM elevation of P (the orthometric height; see the `SRTM collection user guide <https://lpdaac.usgs.gov/sites/default/files/public/measures/docs/NASA_SRTM_V3.pdf>`_) and N is the height of the EGM96 geoid above the WGS84 ellipsoid at P (the geoid height). 
    """
    # Convert to radians and define constants
    lam = radians(lon)
    phi = radians(lat)
    h = height
    lam_s = radians(satellite_lon)
    r = cs.R_S
    a = cs.WGS84_A
    e2 = cs.WGS84_E2
    N = a/sqrt(1 - e2*sin(phi)**2)

    # Transform P and S coordinates from spherical to rectangular
    x_p = (N + h)*cos(lam)*cos(phi)
    y_p = (N + h)*sin(lam)*cos(phi)
    z_p = (N*(1 - e2) + h)*sin(phi)

    x_s = r*cos(lam_s)
    y_s = r*sin(lam_s)
    z_s = 0

    # Translate coordinate system origin to P
    x = x_s - x_p
    y = y_s - y_p
    z = z_s - z_p

    # Transform to P-local geodetic coordinates
    e = -x*sin(lam) + y*cos(lam)
    n = -x*sin(phi)*cos(lam) - y*sin(phi)*sin(lam) + z*cos(phi)
    u = x*cos(phi)*cos(lam) + y*cos(phi)*sin(lam) + z*sin(phi)

    # Compute azimuth and elevation of S relative to P
    alp = atan2(e, n)
    nu = atan2(u, sqrt(e**2 + n**2))

    # Azimuths are positive by convention
    if alp < 0:
        alp += 2*pi

    # Return in degrees
    return degrees(alp), degrees(nu)


def get_geoid_height(lon, lat, num_tries=3):
    """
    Query https://geographiclib.sourceforge.io/cgi-bin/GeoidEval for the height in meters of the EGM96 geoid above the WGS84 ellipsoid for the given longitude and latitude. 
    If the result is negative, then the geoid lies below the ellipsoid.
    Raise a ``ValueError`` if the query fails after ``num_tries`` tries.

    NOTES:
        - It would be good to rewrite this function so that it does not depend on internet access. For starters, see `https://github.com/vandry/geoidheight <https://github.com/vandry/geoidheight>`_, which uses the EGM2008 ellipsoid.
    """
    url = 'https://geographiclib.sourceforge.io/cgi-bin/GeoidEval'
    data = {'input': '{!s}+{!s}'.format(lat, lon)}
    pattern = r'EGM96</a>\s*=\s*<font color="blue">([\d\.\-]+)</font>'

    for i in range(num_tries):
        r = requests.get(url, data)
        if r.status_code != requests.codes.ok:
            continue

        m = re.search(pattern, r.text)
        if m is None:
            raise ValueError('Failed to parse data from', url)
        else:
            return float(m.group(1))

    raise ValueError('Failed to download data from', url)


def compute_satellite_los(in_path, satellite_lon, out_path, n=3,
                          make_shp=False):
    """
    Given the path ``in_path`` to an SRTM1 or SRTM3 file and the longitude of a geostationary satellite, color with 8-bits of grayscale (pixel values from 0 to 255) the raster cells according to whether they are out (blackish, close to 0) or in (whitish, close to 255) of the line-of-site of the satellite, and save the result as a GeoTIFF file located at ``out_path``.
    If ``make_shp``, then also create an ESRI Shapefile bundle (.dbf, .prj, .shp, and .shx files) out of the GeoTIFF and save it to a similar path (same path stem but with Shapefile path suffixes).

    ALGORITHM: 
        #. Partition the SRTM tile into ``n**2`` square subtiles or roughly the same size
        #. For each subtile, compute the longitude, latitude, and (WGS84) height of its center 
        #. Compute the the look angles of the satellite from the center 
        #. Use the look angles to shade the subtile via GDAL's ``gdaldem hillshade`` command
        #. Merge the subtiles and save the result as a GeoTIFF file

    NOTES:
        - Calls :func:`get_geoid_height` ``n**2`` times. Because that function is currently implemented as an HTTP GET request, that slows things down and also introduces ``n**2`` opportunities for failure (raising a ``ValueError``). 
        - To roughly interpret the output raster values as actual satellite signal strengths, one would need to obtain some actual on-the-ground satellite readings.
    """
    in_path = Path(in_path)
    tile_id = ut.get_tile_id(in_path)
    out_path = Path(out_path)
    if not out_path.parent.exists():
        out_path.parent.mkdir(parents=True)

    # Unzip tile file if necessary
    if in_path.name.endswith('.zip'):
        is_zip = True
        shutil.unpack_archive(str(in_path), str(in_path.parent))
        f = in_path.parent/'{!s}.hgt'.format(tile_id)
    else:
        is_zip = False
        f = in_path

    # Create temporary directory to hold the subtiles
    tmp_path = Path(tempfile.mkdtemp())

    # Iterate through subtiles and compute the satellite shadows
    f_info = ut.gdalinfo(f)
    width, height = f_info['width'], f_info['height']
    subtile_names = []
    for i, window in enumerate(partition(width, height, n)):
        # Extract subtile i
        g = tmp_path/'{!s}.tif'.format(i)
        subtile_names.append(g.name)
        args = ['gdal_translate', '-of', 'Gtiff', '-srcwin',
                str(window[0]), str(window[1]), str(window[2]), str(window[3]),
                str(f), str(g)]
        subprocess.run(args, stdout=subprocess.PIPE, universal_newlines=True,
                       check=True)

        # Compute orthometric height H and geoid height N at center of subtile
        lon, lat = ut.gdalinfo(g)['center']
        args = ['gdallocationinfo', str(g), '-wgs84', '-valonly',
                str(lon), str(lat)]
        sp = subprocess.run(args,
                            stdout=subprocess.PIPE, universal_newlines=True, check=True)
        H = float(sp.stdout)
        N = get_geoid_height(lon, lat)

        # Compute look angles at center and then shade with GDAL
        az, el = compute_look_angles(lon, lat, H + N, satellite_lon)
        args = ['gdaldem', 'hillshade', '-compute_edges',
                '-az', str(az), '-alt', str(el), str(g), str(g)]
        subprocess.run(args,
                       stdout=subprocess.PIPE, universal_newlines=True, check=True)

    # Merge subtiles.
    # Use gdalbuildvert and gdal_translate, because gdal_merge.py produces the wrong size image for some reason.
    args = ['gdalbuildvrt', 'merged.vrt'] + [name for name in subtile_names]
    subprocess.run(args, cwd=str(g.parent),
                   stdout=subprocess.PIPE, universal_newlines=True, check=True)

    args = ['gdal_translate', 'merged.vrt', 'merged.tif', '-of', 'GTiff']
    subprocess.run(args, cwd=str(g.parent),
                   stdout=subprocess.PIPE, universal_newlines=True, check=True)

    # Move file
    (g.parent/'merged.tif').replace(out_path)

    # Clean up
    ut.rm_paths(tmp_path)
    if is_zip:
        f.unlink()

    if make_shp:
        tif = out_path.name
        shp = out_path.stem + '.shp'
        args = ['gdal_polygonize.py', str(tif), '-f', 'ESRI Shapefile',
                str(shp)]
        subprocess.run(args, cwd=str(out_path.parent),
                       stdout=subprocess.PIPE, universal_newlines=True, check=True)
