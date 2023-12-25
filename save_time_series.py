import config
from regions_config import wrf_sounding_stations
import wrf_calc
from extra import add_cf_wrf, error_handler, HashDataset
from products import wrf_time_series_products
import xarray as xr
import pandas as pd
import numpy as np
from metpy.units import units
from datetime import datetime
import cartopy.crs as ccrs
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description="Saves WRF time series data.")
parser.add_argument(
    "--date", type=int, required=True,
    help="Starting date of WRF run. (YYYYMMDDHH)")
parser.add_argument(
    "--ens", type=str, default="", help="Ensemble member.")
args = parser.parse_args()

init_time = datetime.strptime(str(args.date), "%Y%m%d%H")


def save_data(data, lon, lat, station_name, init_time, ens):
    point = data["T2"].metpy.cartopy_crs.transform_point(
        lon, lat, ccrs.PlateCarree())
    point_data = data.interp(
        {"west_east": point[0], "south_north": point[1]}, method="linear")
    df = pd.DataFrame({
        "time": [xr.apply_ufunc(
            lambda time: pd.to_datetime(
                time.decode(), format="%Y-%m-%d_%H:%M:%S"),
            point_data["Times"].load(), vectorize=True).dt.strftime(
                "%Y-%m-%d %H:%M").values],
        **dict(zip(
            wrf_time_series_products.keys(),
            [[product.var(point_data).values] for product in
            wrf_time_series_products.values()]))})
    path = f"../images_wrf/{init_time:%Y%m%d%H}/{ens}/time_series/{station_name}/"
    Path(path).mkdir(parents=True, exist_ok=True)
    df.to_json(f"{path}/time_series.json", orient="records")

for domain in np.arange(
        config.wrf_max_domain, config.wrf_min_domain - 1, -1):
    if wrf_sounding_stations[domain] != {}:
        geo = xr.open_dataset(f"../WPS/geo_em.d{domain:02}.nc").squeeze()
        mfds = xr.open_mfdataset(
            f"../wrfout/{init_time:%Y%m%d%H}/{args.ens}/wrfout24_d{domain:02}_*",
            concat_dim="Time", combine="nested")
        mfds = xr.merge(
            [mfds, geo[["COSALPHA", "SINALPHA"]]], compat="override")
        mfds = mfds.drop_dims(["west_east_stag", "south_north_stag"])
        mfds = add_cf_wrf(mfds)
        snapshot = mfds.isel(Time=0)
        snapshot = snapshot.metpy.assign_y_x(tolerance=10 * units("meter"))
        mfds["south_north"], mfds["west_east"] = snapshot["south_north"], snapshot["west_east"]
        mfds["rh2m"] = wrf_calc.rh2m(mfds)
        mfds["umet10"], mfds["vmet10"] = wrf_calc.uvmet10_xr(mfds)
        mfds["umet10"] = mfds["umet10"] * units("m/s")
        mfds["vmet10"] = mfds["vmet10"] * units("m/s")
        for station in wrf_sounding_stations[domain].values():
            save_data(
                mfds, station.lon, station.lat, station.name, init_time,
                args.ens)
