google_bucket = "gs://nbren12-data/data/id/a6243e6261a3532761de9b6eeb14ed3c/10ksamples.nc"

rule raw_data:
    output: "data/raw/rayben/10ksamples.nc"
    shell: "gsutil cp {google_bucket} {output}"


rule regrid_data:
    input: rules.raw_data.output
    output: "data/processed/regridded.nc"
    script: "scripts/regrid.py"

rule stats:
    input: rules.regrid_data.output
    output: mu="data/stats/mu.nc", sig="data/stats/sig.nc"
    script: "scripts/stats.py"

rule bin:
    input: "data/raw/rayben/10ksamples.nc"
    output: "data/processed/binned.nc"
    script: "scripts/bin_isothermal.py"
