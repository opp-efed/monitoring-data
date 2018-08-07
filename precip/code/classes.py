import os
import numpy as np
import pandas as pd

from collections import OrderedDict
from dbfread import DBF, FieldParser


class MemoryMatrix(object):
    """ A wrapper for NumPy 'memmap' functionality which allows the storage and recall of arrays from disk """

    def __init__(self, index, y_size, z_size=None, path=None, base=None, name=None, dtype=np.float32, overwrite=True,
                 input_only=False):
        self.name = name if name else "temp"
        self.dtype = dtype
        self.index = index
        self.count = len(self.index)
        self.lookup = dict(zip(self.index, np.arange(self.count)))
        self.shape = (self.count, y_size) if not z_size else (self.count, y_size, z_size)
        overwrite = False if input_only else overwrite

        # Load from saved file if one is specified, else generate
        default_path = r"C:\Users\Trip Hook\Desktop"  # mkdtemp()
        path = default_path if path is None else path
        if not os.path.exists(path):
            os.makedirs(path)
        base = name if not base else base
        self.path = os.path.join(path, base + ".dat")

        # Determine whether to load existing or create new matrix
        assert not input_only or os.path.exists(self.path), "Matrix {} not found".format(self.path)
        if overwrite or not os.path.exists(self.path):
            np.memmap(self.path, dtype=dtype, mode='w+', shape=self.shape)
            self.new = True
        else:
            self.new = False

    def fetch(self, get_index, verbose=True, raw_index=False, dtype=None, copy=False):
        array = self.copy if copy else self.reader
        location = self.lookup.get(get_index) if not raw_index else get_index
        if location is not None:
            output = array[location]
        else:
            if verbose:
                print("{} not found".format(get_index))
            output = None
        del array
        return output

    def fetch_multiple(self, indices, copy=False, verbose=False, aliased=True, return_index=False, columns=None):

        mode = 'c' if copy else 'r'
        index = None

        if aliased:
            addresses = np.int32([self.lookup.get(x, -1) for x in indices])
            found = np.where(addresses >= 0)[0]
            not_found = indices.size - found.size
            if not_found:
                index = found
                if verbose:
                    print("Missing {} of {} indices in {} matrix".format(not_found, len(addresses), self.name))
            indices = addresses[found]

        array = np.memmap(self.path, dtype='float32', mode=mode, shape=self.shape)
        out_array = array[indices] if columns is None else array[np.ix_(indices, columns)]
        del array

        if return_index:
            return out_array, index
        else:
            return out_array

    def update(self, key, value, aliased=True):
        array = self.writer
        array_index = self.lookup.get(key) if aliased else key
        if array_index is not None:
            array[array_index] = value
        else:
            print("Index {} not found in {} array".format(key, self.name))
        del array

    @property
    def reader(self):
        return np.memmap(self.path, dtype=self.dtype, mode='r', shape=self.shape)

    @property
    def copy(self):
        return np.memmap(self.path, dtype=self.dtype, mode='c', shape=self.shape)

    @property
    def writer(self):
        mode = 'r+' if os.path.isfile(self.path) else 'w+'
        return np.memmap(self.path, dtype=self.dtype, mode=mode, shape=self.shape)


class MetfileMatrix(MemoryMatrix):
    def __init__(self, memmap_path):
        self.dir = os.path.dirname(memmap_path)
        self.name = os.path.basename(memmap_path)
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets
        self.start_date, self.end_date, self.metfiles = self.load_key()

        self.n_dates = int((self.end_date.astype(int) - self.start_date.astype(int))) + 1

        # Set dates
        self.new_years = np.arange(self.start_date, self.end_date + np.timedelta64(365, 'D'),
                                   np.timedelta64(1, 'Y'), dtype='datetime64[Y]').astype('datetime64[D]')
        # Initialize memory matrix
        super(MetfileMatrix, self).__init__(self.metfiles, self.n_dates, 3, path=self.dir,
                                            base=self.name, input_only=True)
    @property
    def dates(self):
        return pd.date_range(self.start_date, self.end_date)

    def load_key(self):

        try:
            data = np.load(self.keyfile_path)
            start_date, end_date = map(np.datetime64, data[:2])
            metfiles = data[2:]
            return start_date, end_date, metfiles
        except ValueError as e:
            print(e)
            exit("Invalid key file {}".format(self.keyfile_path))

    def fetch_station(self, station_id):
        try:
            data = np.array(super(MetfileMatrix, self).fetch(station_id, copy=True, verbose=False)).T
            data[:2] /= 100.  # Precip, PET  cm -> m
            return data
        except:
            print("Met station {} not found".format(station_id))


class Navigator(object):
    def __init__(self, region_id, upstream_path):
        self.file = os.path.join(upstream_path, "region_{}_nav.npz".format(region_id))
        self.paths, self.times, self.map, self.alias_to_reach, self.reach_to_alias = self.load()
        self.reach_ids = set(self.reach_to_alias.keys())

    def load(self):
        assert os.path.isfile(self.file), "Upstream file {} not found".format(self.file)
        data = np.load(self.file, mmap_mode='r')
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion

    def upstream_watershed(self, reach_id, mode='reach', return_times=True):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.reach_to_alias.get(reach_id)
        reaches, adjusted_times, warning = np.array([]), np.array([]), None
        try:
            start_row, end_row, col = map(int, self.map[reach])
            start_col = list(self.paths[start_row]).index(reach)
        except TypeError:
            warning = "Reach {} not found in region".format(reach)
        except ValueError:
            warning = "{} not in upstream lookup".format(reach)
        else:
            # Fetch upstream reaches and times
            aliases = unpack(self.paths)
            reaches = aliases if mode == 'alias' else np.int32(self.alias_to_reach[aliases])
        if not return_times:
            return reaches, warning
        else:
            if warning is None:
                times = unpack(self.times)
                adjusted_times = np.int32(times - self.times[start_row][start_col])
            return reaches, adjusted_times, warning


def read_dbf(dbf_file):
    class MyFieldParser(FieldParser):
        def parse(self, field, data):
            try:
                return FieldParser.parse(self, field, data)
            except ValueError:
                return None

    try:
        dbf = DBF(dbf_file)
        table = pd.DataFrame(iter(dbf))
    except ValueError:
        dbf = DBF(dbf_file, parserclass=MyFieldParser)
        table = pd.DataFrame(iter(dbf))
    table.rename(columns={column: column.lower() for column in table.columns}, inplace=True)
    return table


nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                          ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                          ('03N', {"VA", "NC", "SC", "GA"}),
                          ('03S', {"FL", "GA"}),
                          ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                          ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                          ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                          ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                          ('07', {"MN", "WI", "SD", "IA", "IL", "MO"}),
                          ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                          ('09', {"ND", "MN", "SD"}),
                          ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                          ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                          ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                          ('12', {"NM", "TX", "LA"}),
                          ('13', {"CO", "NM", "TX"}),
                          ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                          ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                          ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                          ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                          ('18', {"OR", "NV", "CA"})))

