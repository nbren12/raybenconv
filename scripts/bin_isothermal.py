"""Bin the data
"""
import dask.array as da
import numpy as np

import xarray as xr


def xcut(x, bins=None):

    # Dask array saves memory
    xd = da.from_array(x, chunks=[100] * len(x.shape))

    # this datatype is good for 2^16 integers
    dtype = np.uint16
    cut = 0 * x.astype(dtype)

    if bins is None:
        bins = np.linspace(0, 1, 20)[1:]

    for b in bins:
        cut.values += np.asarray(x > b)

    return cut


def groupby(x, cuts, f=lambda x: x, post=lambda x: x):
    for c in np.unique(cuts.values):
        print("Cut {cut}".format(cut=c))
        yield post(f(x) * (cuts == c))


def nc2bins(nc, binkey='T', nbins=21, key='uY'):
    """Python binning code

    Parameters
    ----------
    nc : str
        Name of netcdf file with input data. Must contain variables named w, binkey, and key.

    """
    nc = xr.open_dataset(nc)

    print("Generating Cuts")
    bins = np.linspace(0, 1, 20)[1:]
    cuts = xcut(nc[binkey], bins=bins)

    # Average accross bins. If only there was a "groupby" in xr
    print("Averaging accross bins")

    temp = xr.IndexVariable('temp', bins - (bins[1] - bins[0]) / 2)

    def cutmean(nc):
        return nc[key]

    summer = lambda x: x.sum('x')

    binned = xr.concat(groupby(nc, cuts, cutmean, summer), temp)
    return xr.Dataset({key: binned})

    # pop = xr.concat(groupby(nc, cuts, popsum, summer), temp)
    # return xr.Dataset({'w': binned, 'pop': pop})


def main():

    out = nc2bins(snakemake.input[0], nbins=40,
                  key='uY', binkey='T')\
                 .transpose('t','temp','z')

    out['demean'] = out['uY'] - out['uY'].mean('t')

    out.to_netcdf(snakemake.output[0])


try:
    snakemake
except NameError:
    pass
else:
    main()
