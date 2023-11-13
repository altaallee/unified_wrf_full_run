from datetime import datetime, timedelta
import argparse


parser = argparse.ArgumentParser(description="Update namelist files.")
parser.add_argument(
    "--date", type=int, required=True,
    help="Starting date of WRF run. (YYYYMMDDHH)")
parser.add_argument(
    "--hours", type=int, required=True, help="Total hours of run.")
parser.add_argument(
    "--ens", type=str, default=False, help="Ensemble name.")
args = parser.parse_args()

start_date = datetime.strptime(str(args.date), "%Y%m%d%H")
end_date = start_date + timedelta(hours=args.hours)

with open("../WPS/namelist.wps", "r") as file:
    lines = file.readlines()

    for i, line in enumerate(lines):
        line = line.split()
        if line != []:
            if line[0] == "start_date":
                lines[i] = f" start_date = '{start_date:%Y-%m-%d_%H}:00:00','{start_date:%Y-%m-%d_%H}:00:00','{start_date:%Y-%m-%d_%H}:00:00',\n"
            if line[0] == "end_date":
                lines[i] = f" end_date = '{end_date:%Y-%m-%d_%H}:00:00','{start_date:%Y-%m-%d_%H}:00:00','{start_date:%Y-%m-%d_%H}:00:00',\n"

with open("../WPS/namelist.wps", "w") as file:
    file.writelines(lines)

with open(f"../WRF/run{f'_{args.ens}' if args.ens else ''}/namelist.input", "r") as file:
    lines = file.readlines()

    for i, line in enumerate(lines):
        line = line.split()
        if line != []:
            if line[0] == "run_hours":
                lines[i] = f" run_hours = {args.hours},\n"
            if line[0] == "start_year":
                lines[i] = f" start_year = {start_date:%Y}, {start_date:%Y}, {start_date:%Y},\n"
            if line[0] == "start_month":
                lines[i] = f" start_month = {start_date:%m}, {start_date:%m}, {start_date:%m},\n"
            if line[0] == "start_day":
                lines[i] = f" start_day = {start_date:%d}, {start_date:%d}, {start_date:%d},\n"
            if line[0] == "start_hour":
                lines[i] = f" start_hour = {start_date:%H}, {start_date:%H}, {start_date:%H},\n"
            if line[0] == "end_year":
                lines[i] = f" end_year = {end_date:%Y}, {end_date:%Y}, {end_date:%Y},\n"
            if line[0] == "end_month":
                lines[i] = f" end_month = {end_date:%m}, {end_date:%m}, {end_date:%m},\n"
            if line[0] == "end_day":
                lines[i] = f" end_day = {end_date:%d}, {end_date:%d}, {end_date:%d},\n"
            if line[0] == "end_hour":
                lines[i] = f" end_hour = {end_date:%H}, {end_date:%H}, {end_date:%H},\n"

            if line[0] == "history_outname":
                lines[i] = f" history_outname = '../../wrfout/{start_date:%Y%m%d%H}/{args.ens if args.ens else ''}/wrfout_d<domain>_<date>'\n"
            if line[0] == "history_begin":
                lines[i] = f" history_begin = {args.hours * 60 + 1}, {args.hours * 60 + 1}, {args.hours * 60 + 1},\n"
            if line[0] == "history_interval":
                lines[i] = f" history_interval = {args.hours * 60 + 1}, {args.hours * 60 + 1}, {args.hours * 60 + 1},\n"

            if line[0] == "restart_interval":
                lines[i] = f" restart_interval = {args.hours * 60 + 1},\n"

            if line[0] == "auxhist24_outname":
                lines[i] = f" auxhist24_outname = '../../wrfout/{start_date:%Y%m%d%H}/{args.ens if args.ens else ''}/wrfout24_d<domain>_<date>'\n"

with open(f"../WRF/run{f'_{args.ens}' if args.ens else ''}/namelist.input", "w") as file:
    file.writelines(lines)
