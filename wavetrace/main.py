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
import numpy as np
from osgeo import gdal, osr
import wavetrace.constants as cs
import wavetrace.utilities as ut
import math

# In the calls to the subprocess function below,
# sometimes instead of using absolute paths,
# i use relative paths in combination with the ``cwd`` option.
# Mind the difference!


def dbuvm_to_dbm(e_dbuvm: float, f_mhz: float, g_dbi: float = 0.0) -> float:
    """
    60 dBµV/m → ? dBm for a given frequency.

    Formula (ITU-R, FCC, SPLAT! docs):
        P(dBm) = E(dBµV/m) - 77.2 - 20·log10(f_MHz) + G(dBi)
    """
    return e_dbuvm - 77.2 - 20.0 * math.log10(f_mhz) + g_dbi


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
            t['name'] = build_transmitter_name(t['network_name'],
                                               t['site_name'])
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
    Return the text (string) content of a SPLAT! azimuth file (AZ file) corresponding to the given transmitter.

    INPUT:
        - ``transmitter``: dictionary of the same form as any one of the elements in the list output by :func:`read_transmitters`

    OUTPUT:
        String

    NOTES:
        A transmitter with no ``'bearing'`` or ``'horizontal_beamwidth'`` data will produce the string ``'0  0'``.
    """
    t = transmitter
    try:
        bearing = float(t['bearing'])
        hb = float(t['horizontal_beamwidth'])
        left = int(round(360 - (hb/2)))
        right = int(round(hb/2))
        s = '{!s}\n'.format(bearing)
        for x in range(360):
            if left <= x or x <= right:
                normal = 0.9
            else:
                normal = 0.1
            s += '{!s}  {!s}\n'.format(x, normal)
    except:
        s = '0  0\n'

    return s[:-1]  # Drop the final new line


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
                       receiver_sensitivity=cs.RECEIVER_SENSITIVITY,
                       coverage_radius=cs.COVERAGE_RADIUS,
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
    - ``'<transmitter name>.ano'``: Alphanumeric file containing actual signal strength values in dBm

    INPUT:
        - ``in_path``: string or Path object specifying a directory; all the SPLAT! transmitter and elevation data should lie here
        - ``out_path``: string or Path object specifying a directory
        - ``transmitters``: list of transmitter dictionaries (in the form output by :func:`read_transmitters`) to grab the frequencies of to convert dBμV/m to dBm (SPLAT! uses dBm)
        - ``receiver_sensitivity``: float; desired path loss threshold beyond which signal strength contours will not be plotted (measured in dBμV/m)
        - ``coverage_radius``: float; maximum coverage radius in kilometers for SPLAT calculations
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

    # Quick lookup table for transmitter names
    tx_index = {t['name']: t for t in transmitters or []}

    # Get transmitter names
    if transmitters is not None:
        transmitter_names = [t['name'] for t in transmitters]
    else:
        transmitter_names = [p.stem for p in in_path.iterdir()
                             if p.name.endswith('.qth')]

    # Splatify
    splat = 'splat'
    if high_definition:
        splat += '-hd'

    for t in tqdm(transmitter_names, total=len(transmitter_names), desc="Computing coverage"):
        f_mhz = tx_index[t]['frequency']          # each tx has its own freq
        rx_thresh = dbuvm_to_dbm(receiver_sensitivity,
                                 f_mhz)     # 60 dBµV/m → dBm

        print(f"Transmitter {t} has sensitivity {rx_thresh} dBm")

        # build splat argument list - add -ano for alphanumeric output with actual signal values
        args = [splat, '-t', t + '.qth', '-L', '8.0', '-R', str(coverage_radius), '-dbm', '-db',
                str(rx_thresh), '-metric', '-ngs', '-kml', '-ppm', '-ano', t + '.ano',
                '-o', t + '.ppm']
        subprocess.run(args, cwd=str(in_path),
                       stdout=subprocess.PIPE, universal_newlines=True, check=True)

    # Move outputs to out_path
    exts = ['.ppm', '-ck.ppm', '-site_report.txt', '.kml', '.ano']
    for t in transmitter_names:
        for ext in exts:
            src = in_path/(t + ext)
            tgt = out_path/(t + ext)
            shutil.move(str(src), str(tgt))


def create_binary_coverage_from_ano(ano_file, kml_file, threshold_dbm, output_tiff):
    """
    Read a SPLAT! alphanumeric output file (.ano) containing signal strength values in dBm,
    apply a precise threshold, and create a binary GeoTIFF (0/1 values) where 1 indicates
    signal strength >= threshold_dbm.

    INPUT:
        - ``ano_file``: Path to SPLAT! .ano file containing actual dBm signal values
        - ``kml_file``: Path to corresponding .kml file for geographic bounds
        - ``threshold_dbm``: float; signal strength threshold in dBm
        - ``output_tiff``: Path where binary GeoTIFF will be saved

    OUTPUT:
        None. Creates a binary GeoTIFF file.

    NOTES:
        - SPLAT! .ano files have format: lat, lon, azimuth, elevation, signal_strength_dbm
        - Lines ending with '*' indicate obstructed paths
        - Only areas with signal >= threshold_dbm will have value 1, others 0
    """
    ano_file = Path(ano_file)
    kml_file = Path(kml_file)
    output_tiff = Path(output_tiff)

    # Read geographic bounds from KML
    with kml_file.open() as f:
        kml_content = f.read()
    bounds = get_bounds_from_kml(kml_content)
    min_lon, min_lat, max_lon, max_lat = bounds

    # Parse .ano file
    lats, lons, signals = [], [], []

    with ano_file.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith(';') or ',' not in line:
                continue  # Skip comments and header lines

            parts = line.rstrip('*').split(',')  # Remove '*' marker and split
            if len(parts) >= 5:
                try:
                    lat = float(parts[0])
                    lon = float(parts[1])
                    # 5th column is signal strength
                    signal_dbm = float(parts[4])

                    lats.append(lat)
                    lons.append(lon)
                    signals.append(signal_dbm)
                except ValueError:
                    continue  # Skip malformed lines

    if not lats:
        raise ValueError(f"No valid data found in {ano_file}")

    # Convert to numpy arrays
    lats = np.array(lats)
    lons = np.array(lons)
    signals = np.array(signals)

    # Apply threshold: 1 for signal >= threshold_dbm, 0 otherwise
    binary_values = (signals >= threshold_dbm).astype(np.uint8)

    # Create raster grid
    # Determine grid resolution based on data density
    lat_range = max_lat - min_lat
    lon_range = max_lon - min_lon

    # Use reasonable resolution - adjust if needed
    # ~1 arc-second resolution
    n_lat = max(100, min(1000, int(lat_range * 3600)))
    n_lon = max(100, min(1000, int(lon_range * 3600)))

    lat_step = lat_range / n_lat
    lon_step = lon_range / n_lon

    # Create grid
    raster = np.zeros((n_lat, n_lon), dtype=np.uint8)

    # Fill raster with binary values
    for lat, lon, val in zip(lats, lons, binary_values):
        if min_lat <= lat <= max_lat and min_lon <= lon <= max_lon:
            row = int((max_lat - lat) / lat_step)
            col = int((lon - min_lon) / lon_step)

            # Ensure indices are within bounds
            row = max(0, min(n_lat - 1, row))
            col = max(0, min(n_lon - 1, col))

            raster[row, col] = val

    # Create GeoTIFF
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(str(output_tiff), n_lon, n_lat, 1, gdal.GDT_Byte)

    # Set geotransform (affine transformation)
    geotransform = [min_lon, lon_step, 0, max_lat, 0, -lat_step]
    dataset.SetGeoTransform(geotransform)

    # Set projection (WGS84)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    dataset.SetProjection(srs.ExportToWkt())

    # Write data
    band = dataset.GetRasterBand(1)
    band.WriteArray(raster)
    band.SetNoDataValue(0)  # Set 0 as no-data value

    # Clean up
    dataset = None


def postprocess_coverage_0(path, keep_ppm, make_shp):
    """
    Using the ANO and KML files in the directory ``path`` do the following:

    - Create binary GeoTIFF files from ANO files using precise signal thresholds
    - Convert binary GeoTIFF to PNG for visualization
    - Update KML files to reference PNG files
    - Optionally create ESRI Shapefile bundles (.dbf, .prj, .shp, and .shx files) from the GeoTIFF files

    INPUT:
        - ``path``: string or Path object; directory where coverage reports (outputs of :func:`compute_coverage`) lie
        - ``keep_ppm``: boolean; keep the original, large PPM files in the coverage reports if and only if this flag is ``True``
        - ``make_shp``: boolean; create ESRI Shapefiles from the GeoTIFF files if and only if this flag is ``True``

    OUTPUT:
        None.

    NOTES:
        - This new approach reads the .ano files containing actual dBm values
        - Creates precise binary thresholds instead of relying on color-coded PPM files
        - Assumes the threshold was already applied during SPLAT! execution via the -db flag
    """
    path = Path(path)

    # First pass: create binary GeoTIFF from ANO files
    ano_files = list(path.glob('*.ano'))
    for ano_file in tqdm(ano_files, desc="Creating binary coverage maps"):
        transmitter_name = ano_file.stem
        kml_file = path / f"{transmitter_name}.kml"
        tiff_file = path / f"{transmitter_name}.tif"

        if kml_file.exists():
            # The threshold was already applied by SPLAT! via the -db flag
            # So any signal in the .ano file is >= our threshold
            # We'll use a very low threshold (-200 dBm) to capture all reported signals
            create_binary_coverage_from_ano(
                ano_file, kml_file, -200.0, tiff_file)

    # Second pass: create PNG from GeoTIFF for visualization
    for tiff_file in path.glob('*.tif'):
        if tiff_file.name.endswith('-ck.tif'):
            continue  # Skip legend files

        png_file = tiff_file.with_suffix('.png')

        # Convert binary GeoTIFF to PNG
        # Use white for coverage (value 1) and transparent for no coverage (value 0)
        args = ['gdal_translate', '-of', 'PNG', '-scale', '0', '1', '0', '255',
                str(tiff_file), str(png_file)]
        subprocess.run(args, cwd=str(path),
                       stdout=subprocess.PIPE, universal_newlines=True, check=True)

        # Make value 0 (no coverage) transparent
        args = ['convert', '-transparent',
                '#000000', str(png_file), str(png_file)]
        subprocess.run(args, cwd=str(path),
                       stdout=subprocess.PIPE, universal_newlines=True, check=True)

    # Third pass: update KML files to reference PNG instead of PPM
    for kml_file in path.glob('*.kml'):
        with kml_file.open() as f:
            kml_content = f.read()

        # Replace .ppm with .png in KML
        kml_content = kml_content.replace('.ppm', '.png')

        with kml_file.open('w') as f:
            f.write(kml_content)

    # Fourth pass: clean up PPM files if requested
    if not keep_ppm:
        for ppm_file in path.glob('*.ppm'):
            ppm_file.unlink()

    # Optional fifth pass: create vector files from GeoTIFFs
    if make_shp:
        for tiff_file in path.glob('*.tif'):
            if tiff_file.name.endswith('-ck.tif'):
                continue  # Skip legend files

            shp_file = tiff_file.with_suffix('.shp')
            args = ['gdal_polygonize.py', str(
                tiff_file), '-f', 'ESRI Shapefile', str(shp_file)]
            subprocess.run(args, cwd=str(path),
                           stdout=subprocess.PIPE, universal_newlines=True, check=True)


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
                     coverage_radius=cs.COVERAGE_RADIUS,
                     keep_ppm=False, high_definition=False, make_shp=False):
    """
    Produce coverage reports by running :func:`compute_coverage_0` and then run :func:`post_process_coverage_0`.
    """
    compute_coverage_0(in_path, out_path, transmitters, receiver_sensitivity,
                       coverage_radius, high_definition)
    postprocess_coverage_0(out_path, keep_ppm, make_shp)


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
