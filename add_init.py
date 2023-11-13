import pandas as pd
from datetime import datetime
import argparse


parser = argparse.ArgumentParser(description="Add init time to json file.")
parser.add_argument(
    "--date", type=int, help="Date of init time formatted as YYYYMMDDHH.")
parser.add_argument(
    "--hours", type=int, help="Total hours of run.")
parser.add_argument(
    "--ens", type=str, default="", help="Ensemble name.")
args = parser.parse_args()

path = "init_times.json"

df = pd.read_json(path)
df = pd.concat(
    [df, pd.DataFrame(
        {"name": [f"""{datetime.strftime(
            datetime.strptime(str(args.date), '%Y%m%d%H'), '%HZ %b %d')}{f' - {args.ens}' if args.ens != '' else ''}"""], 
         "value": [args.date], "fcsthours": [args.hours], "ens": [args.ens]})],
        axis=0)
df.sort_values(["value", "ens"], ascending=False, inplace=True)
df.to_json(path, orient="records")
