import numpy as np
from scipy.interpolate import interp2d
import xarray as xr
import toolz


def get_regular_grid(nx, nz):
    xg = np.linspace(0, 15.0, nx, endpoint=False)
    zg = np.linspace(0, 1.0, nz, endpoint=False)

    return xg, zg


@toolz.curry
def interpolate_data(x, z, xg, zg, data):
    return interp2d(x, z, data[0, :, :])(xg, zg)[np.newaxis, :, :]


def remove_duplicate_x(ds):
    x = np.asarray(ds.x)
    xind, = np.nonzero(np.diff(x) > 1e-10)
    return ds.isel(x=xind)


def regrid_dataarray(xarr):

    ds_deduped = remove_duplicate_x(xarr)
    darr = ds_deduped.data

    x = np.asarray(ds_deduped.x)
    z = np.asarray(ds_deduped.z)

    xg, zg = get_regular_grid(512, 64)

    interpolater = interpolate_data(x, z, xg, zg)
    # map_blocks
    darr_interp = darr.map_blocks(
        interpolater, dtype=darr.dtype, chunks=(1, len(zg), len(xg)))
    out_xarr = xr.DataArray(
        darr_interp,
        coords={'x': xg,
                'z': zg,
                't': ds_deduped['t']},
        dims=['t', 'z', 'x'])

    return out_xarr


def regrid_dataset(ds):
    return ds.apply(regrid_dataarray)


xr.open_dataset(snakemake.input[0], chunks={'t': 1})\
  .apply(regrid_dataarray)\
  .to_netcdf(snakemake.output[0])
