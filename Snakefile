google_bucket = "gs://nbren12-data/data/id/a6243e6261a3532761de9b6eeb14ed3c/10ksamples.nc"

rule raw_data:
    output: "data/raw/rayben/10ksamples.nc"
    shell: "gsutil cp {google_bucket} {output}"
