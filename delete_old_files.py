import pandas as pd
from datetime import datetime, timedelta
import subprocess
import argparse


parser = argparse.ArgumentParser(
    description="Delete files over a given number of days old.")
parser.add_argument(
    "--days", type=int, default=14, 
    help="Number of days old files will be deleted.")
args = parser.parse_args()

path = "init_times.json"
df = pd.read_json(path)

for i, row in df.iterrows():
    init_date = datetime.strptime(str(row.value), "%Y%m%d%H")
    if (datetime.now() - init_date > timedelta(days=args.days)):
        df.drop(df.index[i], inplace=True)
        subprocess.run(f"rm -rf ../images_wrf/{init_date:%Y%m%d%H}", shell=True)
        subprocess.run(f"rm -rf ../wrfout/{init_date:%Y%m%d%H}", shell=True)

df.sort_values(["value"], ascending=False, inplace=True)
df.to_json(path, orient="records")
