# loading the DDOR 2014 Radio Source Catalog
#
#  Trying to keep this similar to the skyfield hipparcos.py load procedure...

URL = 'https://deepspace.jpl.nasa.gov/dsndocs/810-005/107/catalog-fixed.txt'

url = URL

PANDAS_MESSAGE = """Skyfield needs Pandas to load the DDOR 2014 Radio Source Catalog

To load this catalog, Skyfield needs the Pandas data
analysis toolkit. Try installing it using your usual Python package
installer, like "pip install pandas" or "conda install pandas".
"""

_COLUMN_NAMES = (
    'B1950', 'Common', 'No.', 'H', 'M', 'S', 'n',
    'D', 'M', 'S', '(s)', '(asec)',
    'Corr.', 'First', 'Last', 'Obs.',
    '1060', '1040', 'Ind.',
    )

def load_dataframe(fobj):
    """Given an open file for 'catalog-fixed.txt', return a parsed dataframe.

    """
    try:
        from pandas import read_fwf
    except ImportError:
        raise ImportError(PANDAS_MESSAGE)

    fobj.seek(0)
    magic = fobj.read(2)
    compression = 'gzip' if (magic == b'\x1f\x8b') else None
    fobj.seek(0)

    df = read_fwf(
            fobj, compression=compression,
            skiprows=22,
            colspecs=[(1,9),(11,22),(24,28),(30,31),(34,35),(38,49),(50,51),
                        (52,53),(56,57),(59,68),(122,125),(128,131),
                        ]
            )
    df.columns = (
            'B1950', 'name', 'id', 'ra_hrs', 'ra_minutes', 'ra_seconds',
            'sign', 'dec_deg', 'dec_minutes', 'dec_seconds', 'flux_1060',
            'flux_1040',
            )
    df = df.assign(
            ra_hrs = df['ra_hrs'].fillna(0),
            ra_minutes = df['ra_minutes'].fillna(0),
            dec_minutes = df['dec_minutes'].fillna(0),
            dec_seconds = df['dec_seconds'].fillna(0),
            epoch_year = 2000.0,
            )
    df['ra_hours'] = df.apply(ra_to_set, axis=1)
    df['dec_degrees'] = df.apply(dec_to_set, axis=1)

    return df.set_index('name')

def dec_to_set(row):
    """Convert the dec columns into a set, and fix the sign.
    """
    if '-' in row['sign']:
        return (-1*row['dec_deg'],row['dec_minutes'],row['dec_seconds'])
    else:
        return (row['dec_deg'],row['dec_minutes'],row['dec_seconds'])

def ra_to_set(row):
    """Convert the ra columns into a set.
    """
    return (row['ra_hrs'],row['ra_minutes'],row['ra_seconds'])


