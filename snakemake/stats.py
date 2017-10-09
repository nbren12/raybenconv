import xarray as xr

ds = xr.open_dataset(snakemake.input[0], chunks={'t':1000})
mu = ds.mean(['t', 'x'])
sig = ds.std(['t', 'x'])
mu.to_netcdf(snakemake.output.mu)
sig.to_netcdf(snakemake.output.sig)
