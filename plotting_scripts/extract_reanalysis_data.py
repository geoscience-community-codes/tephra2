import argparse
import xarray
import re
import numpy


triple = "([0-9]+)(?:.)([0-9]+)(?:')([0-9\.]+)(?:\")"
PARSER = {
    'lat': re.compile(triple+"([NS])"),
    'lon': re.compile(triple+"([EW])")}


def parse_location(loc, axis='lat'):
    try:
        location = float(loc)
    except ValueError:
        m = PARSER[axis].match(loc)
        if m is None:
            raise RuntimeError(f'could not parser {axis} {loc}')
        location = float(m[1]) \
            + float(m[2]) / 60 \
            + float(m[3]) / 3600
        if m[4] in ['S', 'W']:
            location = -location
    return location


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', help="input netCDF files")
    parser.add_argument('-d', '--date', required=True,
                        help="the date for which to extract data")
    parser.add_argument('-l', '--location', nargs=2, metavar=("LAT", "LON"),
                        required=True,
                        help="the location for which to extract data")
    parser.add_argument('-o', '--output', default="wind",
                        help="name of output file, default=wind")
    args = parser.parse_args()

    dataset = xarray.open_mfdataset(args.input, combine='by_coords')

    ok = True
    for v in ['hgt', 'uwnd', 'vwnd']:
        if v not in dataset.variables:
            ok = False
            print(f'variable {v} not in dataset')
    if not ok:
        parser.error('not all required variables in dataset')

    try:
        lat = parse_location(args.location[0], axis='lat')
        lon = parse_location(args.location[1], axis='lon')
    except RuntimeError as e:
        parser.error(e)
    if lon < 0:
        lon = 180+lon

    print(lat, lon, args.date)

    # select column of data
    col = dataset.sel(time=args.date, lat=lat, lon=lon, method='nearest')
    col = col.mean(dim="time")

    # extract wind speed and direction
    speed = (col.uwnd * col.uwnd + col.vwnd * col.vwnd)**0.5
    direction = numpy.degrees(numpy.arctan2(col.vwnd, col.uwnd))

    with open(args.output, 'w') as out:
        out.write("#HEIGHT	SPEED	DIRECTION\n")
        for l in range(col.level.size):
            out.write(
                f'{col.hgt[l].values} {speed[l].values} '
                f'{direction[l].values}\n')


if __name__ == '__main__':
    main()
