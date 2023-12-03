import pandas as pd
import multiprocessing as mp
import config
from regions_config import wrf_sounding_stations
import argparse
from datetime import datetime, timedelta
import wrf_calc
from extra import add_cf_wrf, error_handler, read_data_wrf_sounding
import cartopy.crs as ccrs
from metpy.units import units
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from pathlib import Path
import numpy as np
import matplotlib.patheffects as path_effects
import matplotlib.transforms as transforms


parser = argparse.ArgumentParser(description="Plot skewt of WRF output.")
parser.add_argument(
    "--date", type=int, required=True,
    help="Starting date of WRF run. (YYYYMMDDHH)")
parser.add_argument(
    "--hours", type=int, required=True, help="Total hours to plot data.")
parser.add_argument(
    "--ens", type=str, default="", help="Ensemble member.")
args = parser.parse_args()

start_date = datetime.strptime(str(args.date), "%Y%m%d%H")
end_date = start_date + timedelta(hours=args.hours)


class Skewt():
    def __init__(self, figsize=[12, 8]):
        mosaic = [
            ["omega", "skewt", "hodo", "theta"],
            ["omega", "skewt", "indices", "indices"],
        ]
        self.fig, self.axes = plt.subplot_mosaic(
            mosaic, width_ratios=[0.75, 7, 2.7, 2], height_ratios=[3, 4],
            gridspec_kw={
                "hspace": 0, "wspace": 0.2, "left": 0.025, "right": 0.975},
            figsize=figsize)

        # skewt
        ss = self.axes["skewt"].get_subplotspec()
        self.axes["skewt"].remove()
        self.axes["skewt"] = self.fig.add_subplot(ss, projection="skewx")
        self.skewt = SkewT(
            fig=self.fig, subplot=self.axes["skewt"], aspect=100)
        self.skewt.plot_dry_adiabats(linewidth=0.5)
        self.skewt.plot_mixing_lines(linewidth=0.5)
        self.skewt.plot_moist_adiabats(linewidth=0.5)
        self.skewt.ax.axvline(-20, color="blue", linewidth=0.5, linestyle="--")
        self.skewt.ax.axvline(0, color="blue", linewidth=0.5, linestyle="--")

        # omega
        skewt_pos = self.skewt.ax.get_position()
        omega_pos = self.axes["omega"].get_position()
        self.axes["omega"].set_position(
            (omega_pos.x0, skewt_pos.y0, omega_pos.width, skewt_pos.height))
        self.axes["omega"].set_yscale("log")
        self.axes["omega"].sharey(self.skewt.ax)
        self.axes["omega"].set_xlim(1, -6)
        self.axes["omega"].set_xticks([0, -2, -4])
        self.axes["omega"].set_xlabel("Omega\nPa s$^{-1}$")
        self.axes["omega"].get_yaxis().set_visible(False)
        self.axes["omega"].spines["top"].set_visible(False)
        self.axes["omega"].spines["right"].set_visible(False)
        self.axes["omega"].spines["bottom"].set_visible(False)
        self.axes["omega"].spines["left"].set_visible(False)

        # hodo
        hodo_pos = self.axes["hodo"].get_position()
        theta_pos = self.axes["theta"].get_position()
        self.axes["hodo"].set_position(
            (hodo_pos.x0, hodo_pos.y0, theta_pos.x0 - hodo_pos.x0,
             hodo_pos.height))
        self.hodo = Hodograph(self.axes["hodo"])
        self.hodo.add_grid(increment=20, linestyle="-", linewidth=0.5)
        for i in [20, 40, 60, 80]:
            self.hodo.ax.text(
                i, 0, f"\n{i}", ha="center", va="center", color="grey",
                fontsize=6)
        self.hodo.ax.set_xticks([])
        self.hodo.ax.set_yticks([])

        # theta
        self.axes["theta"].tick_params(axis="y",direction="in", pad=-25)
        self.axes["theta"].tick_params(axis="x",direction="in", pad=-15)
        self.axes["theta"].set_ylim(999, 501)

        # indices
        self.axes["indices"].set_xlim(0, 1)
        self.axes["indices"].set_ylim(0, 1)
        self.axes["indices"].axis("off")

    def plot_1_6km_barbs(self, u_1km, v_1km, u_6km, v_6km):
        """
        Plots wind barbs at 1 and 6km.

        Parameters
        ----------
        u_1km : float
            U component of wind at 1km in knots.
        v_1km : float
            v component of wind at 1km in knots.
        u_6km : float
            U component of wind at 6km in knots.
        v_6km : float
            v component of wind at 6km in knots.
        """
        self.axes["indices"].barbs(0.85, 0.6, u_1km, v_1km, color="red")
        self.axes["indices"].barbs(0.85, 0.6, u_6km, v_6km, color="blue")

    def plot_barbs(self, pressure, u, v, length=6, linewidth=1, **kwargs):
        """
        Pots wind barbs on side of skewt.

        Parameters
        ----------
        pressure : list
            Pressure values in millibars.
        u : list
            U component of wind in knots.
        v : list
            V component of wind in knots.
        length : int (default=6)
            Length of wind barbs.
        linewidth : int (default=1)
            Thickness of wind barbs.
        """
        p_barb = [i for i in np.arange(1000, 99, -50) if
                  (i > pressure[-1]) & (i < pressure[0])]
        u_barb = wrf_calc.interpolate_1d(np.log(pressure), u, np.log(p_barb))
        v_barb = wrf_calc.interpolate_1d(np.log(pressure), v, np.log(p_barb))
        self.skewt.plot_barbs(
            p_barb, u_barb, v_barb, length=length, linewidth=linewidth,
            **kwargs)

    def plot_critical_angle(self, angle):
        """
        Plots critical angle value in hodograph.

        Parameters
        ----------
        angle : float
            Critical angle.
        """
        self.hodo.ax.text(
            0, 0, f" Critical Angle: {utils.INT2STR(angle)}", ha="left",
            va="bottom", fontsize=8, transform=self.hodo.ax.transAxes)

    def plot_dgz(self, bottom_prs, top_prs):
        """
        Plots DGZ heights on omega plot.

        Parameters
        ----------
        bottom_prs : float
            Pressure at bottom of DGZ in millibars.
        top_prs : float
            Pressure at top of DGZ in millibars.
        """
        self.axes["omega"].axhline(
            bottom_prs, color="cyan", linewidth=3, linestyle="--")
        self.axes["omega"].axhline(
            top_prs, color="cyan", linewidth=3, linestyle="--")

    def plot_heights(self, pressure, h):
        """
        Plots height indicators on skewt.

        Parameters
        ----------
        pressure : list
            Pressure values in millibars.
        h : list
            Height values in meters.
        """
        heights = [i for i in [1, 3, 6, 9, 12, 15] if
                   (i < max(h) / 1000) & (i > min(h) / 1000)]
        prs = wrf_calc.interpolate_1d(h / 1000, pressure, heights)
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.text(
            0.01, pressure.max(), f"SFC", color="red", ha="left", va="center",
            fontsize=8, transform=trans)
        for p, z in zip(prs, heights):
            self.skewt.ax.text(
                0.01, p, f"{z} km", color="red", ha="left", va="center",
                fontsize=8, transform=trans)

    def plot_hodo(self, height, u, v):
        """
        Plots hodograph line.

        Parameters
        ----------
        height : list
            Height values in kilometers.
        u : list
            U component of wind in knots.
        v : list
            V component of wind in knots.
        """
        height_new = height[height < 15]
        u = u[height < 15]
        v = v[height < 15]
        self.hodo.plot_colormapped(
            u, v, height_new * units(""),
            intervals=np.array([0, 3, 6, 9, 99]) * units(""),
            colors=["red", "green", "orange", "blue"], linewidth=1.5)
        heights_label = [
            i for i in [1, 2, 3, 6, 9] if (i > height[0]) & (i < height[-1])]
        u_label = wrf_calc.interpolate_1d(height_new, u, heights_label)
        v_label = wrf_calc.interpolate_1d(height_new, v, heights_label)
        for h, x, y in zip(heights_label, u_label, v_label):
            txt = self.hodo.ax.text(x, y, f"{h}", ha="center", va="center")
            txt.set_path_effects(
                [path_effects.Stroke(linewidth=2, foreground="white"),
                    path_effects.Normal()])

    def plot_hodo_inflow(self, bottom_u, bottom_v, top_u, top_v, movement_u,
                         movement_v):
        """
        Plots inflow layer on hodograph.

        Parameters
        ----------
        bottom_u : float
            U component of wind at bottom of inflow layer in knots.
        bottom_v : float
            V component of wind at bottom of inflow layer in knots.
        top_u : float
            U component of wind at top of inflow layer in knots.
        top_v : float
            V component of wind at top of inflow layer in knots.
        movement_u : float
            U component of storm motion.
        movement_v : float
            V component of storm motion.
        """
        self.hodo.ax.plot(
            [bottom_u, movement_u, top_u], [bottom_v, movement_v, top_v],
            color="cyan")

    def plot_indices(self, indices):
        """
        Plots indices and kinetics.

        Parameters
        ----------
        indices : array
            Array formatted as table with indices and kinetics.
        """
        table = self.axes["indices"].table(
            indices, cellLoc="center", loc="center", edges="")
        table.auto_set_font_size(False)
        table.set_fontsize(9)

    def plot_inflow(self, bottom_prs, top_prs, bottom_hgt, top_hgt, srh):
        """
        Plots inflow layer on skewt.

        Parameters
        ----------
        bottom_prs : float
            Bottom pressure of inflow layer in millibars.
        top_prs : float
            Top pressure of inflow layer in millibars.
        bottom_hgt : float
            Bottom height of inflow layer in meters.
        top_hgt : float
            Top height of inflow layer in meters.
        srh : float
            SRH of inflow layer.
        """
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [0.4, 0.5], [top_prs, top_prs], color="purple", transform=trans)
        self.skewt.ax.plot(
            [0.45, 0.45], [bottom_prs, top_prs], color="purple",
            transform=trans)
        self.skewt.ax.plot(
            [0.4, 0.5], [bottom_prs, bottom_prs], color="purple",
            transform=trans)
        self.skewt.ax.text(
            0.4, top_prs, f"{utils.INT2STR(top_hgt)} m   ", fontsize=8,
            ha="right", va="bottom", color="purple", transform=trans)
        self.skewt.ax.text(
            0.4, bottom_prs, f"{utils.INT2STR(bottom_hgt)} m   ",
            fontsize=8, ha="right", va="top", color="purple", transform=trans)
        self.skewt.ax.text(
            0.45, top_prs, f"{utils.INT2STR(srh)} m$^{2}$ s$^{{{-2}}}$",
            fontsize=8, ha="center", va="bottom", color="purple",
            transform=trans)

    def plot_lines_skewt(self, pressure, t, **kwargs):
        """
        Plots lines on skewt.

        Parameters
        ----------
        pressure : list
            Pressure values in millibars.
        t : list
            Temperature values in celsius.
        """
        self.skewt.plot(pressure, t, **kwargs)

    def plot_max_lr(self, bottom_prs, top_prs, lr):
        """
        Plots maximum lapse rate indicator on skewt.

        Parameters
        ----------
        bottom_prs : float
            Bottom pressure of maximum lapse rate in millibars.
        top_prs : float
            Top pressure of maximum lapse rate in millibars.
        lr : float
            Lapse rate of layer in K/km.
        """
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [0.75, 0.8], [top_prs, top_prs], color="peru", transform=trans)
        self.skewt.ax.plot(
            [0.775, 0.775], [bottom_prs, top_prs], color="peru",
            transform=trans)
        self.skewt.ax.plot(
            [0.75, 0.8], [bottom_prs, bottom_prs], color="peru",
            transform=trans)
        self.skewt.ax.text(
            0.775, top_prs, f"{utils.FLOAT2STR(lr, 1)} °C km$^{{{-1}}}$",
            fontsize=8, ha="center", va="bottom", color="peru", transform=trans)

    def plot_omega(self, pressure, omega, **kwargs):
        """
        Plots omega values at side of skewt.

        Parameters
        ----------
        pressure : list
            Pressure values in millibars.
        omega : list
            Omega values in Pa / s
        """
        width = lambda p, w: 10**(np.log10(p) + w / 2) - 10**(np.log10(p) - w / 2)
        mask_pos = (omega > 0) & (pressure > 100)
        mask_neg = (omega < 0) & (pressure > 100)
        self.axes["omega"].barh(
            pressure[mask_pos], omega[mask_pos], color="blue",
            height=width(pressure[mask_pos], 0.01), clip_on=False, **kwargs)
        self.axes["omega"].barh(
            pressure[mask_neg], omega[mask_neg], color="orange",
            height=width(pressure[mask_neg], 0.01), clip_on=False, **kwargs)

    def plot_sounding_levels(self, lcl, lfc, el):
        """
        Plots LCL, LFC, and EL levels on skewt.

        Parameters
        ----------
        lcl : float
            Pressure of LCL in millibars.
        lfc : float
            Pressure of LFC in millibars.
        el : float
            Pressure of EL in millibars.
        """
        trans = transforms.blended_transform_factory(
            self.skewt.ax.transAxes, self.skewt.ax.transData)
        self.skewt.ax.plot(
            [0.85, 0.9], [lcl, lcl], color="lime", linewidth=2, transform=trans)
        self.skewt.ax.text(
            0.875, lcl, "LCL", color="lime", ha="center", va="top",
            transform=trans)
        self.skewt.ax.plot(
            [0.85, 0.9], [lfc, lfc], color="orange", linewidth=2,
            transform=trans)
        self.skewt.ax.text(
            0.875, lfc, "LFC", color="orange", ha="center", va="bottom",
            transform=trans)
        self.skewt.ax.plot(
            [0.85, 0.9], [el, el], color="purple", linewidth=2, transform=trans)
        self.skewt.ax.text(
            0.875, el, "EL", color="purple", ha="center", va="bottom",
            transform=trans)

    def plot_storm_movers(self, right_u, right_v, left_u, left_v, mean_u,
                          mean_v):
        """
        Plots left, right, and mean storm mortion on hodograph.

        Parameters
        ----------
        right_u : float
            U component of motion for right moving storm in knots.
        right_v : float
            V component of motion for right moving storm in knots.
        left_u : float
            U component of motion for left moving storm in knots.
        left_v : float
            V component of motion for left moving storm in knots.
        mean_u : float
            U component of mean motion for storm in knots.
        mean_v : float
            v component of mean motion for storm in knots.
        """
        self.hodo.ax.plot(
            [right_u, left_u], [right_v, left_v], linestyle="", marker="o",
            fillstyle="none", color="black")
        self.hodo.ax.plot(
            [mean_u], [mean_v], linestyle=None, marker="s", fillstyle="none",
            color="brown")
        self.hodo.ax.text(
            right_u, right_v, f"  RM", fontsize=8, ha="left", va="center")
        self.hodo.ax.text(
            left_u, left_v, f"  LM", fontsize=8, ha="left", va="center")

    def plot_surface_values(self, p, t, td):
        """
        Plots surface values at bottom of skewt.

        Parameters
        ----------
        p : float
            Pressure at surface in millibars.
        t : float
            Temperature at surface in celsius.
        td : float
            Dew point at surface in celsius.
        """
        self.skewt.ax.text(
            td, p, f"{utils.INT2STR(thermo.ctof(td))}°F", ha="center", va="top",
            color="green", weight="bold")
        self.skewt.ax.text(
            t, p, f"{utils.INT2STR(thermo.ctof(t))}°F", ha="center", va="top",
            color="red", weight="bold")

    def plot_thetae(self, pressure, thetae):
        """
        Plots theta profile.

        Parameters
        ----------
        pressure : list
            Pressure values in millibars.
        theta : list
            Theta values in kelvin.
        """
        mask = pressure > 500
        self.axes["theta"].plot(thetae[mask], pressure[mask], color="blue")

    def plot_title(self, station_name, lon, lat, fcst_time, init_time):
        """
        Plots title and time into on skewt.

        Parameters
        ----------
        station_name : str
            Name of station.
        lon : float
            Longitude of station.
        lat : float
            Latitude of station.
        fcst_time : datetime.datetime
            Forecast time of skewt.
        init_time : datetime.datetime
            Initalization time of model.
        """
        self.axes["skewt"].set_title(
            f"{station_name} Lon:{round(lon, 2)} Lat:{round(lat, 2)}",
            loc="left")
        self.axes["theta"].set_title(
            f"Init: {init_time:%Y-%m-%d %H:%MZ}\nValid: {fcst_time:%Y-%m-%d %H:%MZ}",
            loc="right")

    def save_image(self, path, fname, dpi=100, facecolor="white", **kwargs):
        """
        Saves image as PNG.

        Parameters
        ----------
        path : str
            Path of image.
        fname : str
            Name of image.
        dpi : int (default=100)
            The resolution in dots per inch. If 'figure', use the figure's dpi
            value.
        facecolor : str (default="white")
            The facecolor of the figure. If 'auto', use the current figure
            facecolor.
        """
        Path(path).mkdir(parents=True, exist_ok=True)
        self.fig.savefig(
            f"{path}/{fname}.png", dpi=dpi, facecolor=facecolor, **kwargs)
        plt.close()


def plot_sounding(data, lon, lat, station_name, fcst_time, init_time, ens):
    point = data["T2"].metpy.cartopy_crs.transform_point(
        lon, lat, ccrs.PlateCarree())
    data = data.metpy.assign_y_x(tolerance=10 * units("meter"))
    point_data = data.interp(
        {"west_east": point[0], "south_north": point[1]}, method="linear")
    prof = profile.create_profile(
        profile="convective", pres=wrf_calc.pressure(point_data) / 100,
        hght=wrf_calc.height(point_data),
        tmpc=thermo.ktoc(wrf_calc.temperature(point_data)),
        dwpc=wrf_calc.dewpoint(point_data),
        u=utils.MS2KTS(point_data["umet"]),
        v=utils.MS2KTS(point_data["vmet"]),
        omeg=point_data["omega"])

    skewt = Skewt()

    skewt.plot_title(station_name, lon, lat, fcst_time, init_time)

    # omega
    skewt.plot_omega(prof.pres, prof.omeg)
    skewt.plot_dgz(prof.dgz_pbot, prof.dgz_ptop)

    # skewt
    skewt.plot_lines_skewt(prof.pres, prof.wetbulb, color="cyan", linewidth=1)
    skewt.plot_lines_skewt(
        prof.dpcl_ptrace, prof.dpcl_ttrace, color="purple", linewidth=1,
        linestyle="--")
    skewt.plot_lines_skewt(prof.pres, prof.dwpc, color="green", linewidth=2)
    skewt.plot_lines_skewt(
        prof.pres, prof.vtmp, color="red", linewidth=1, linestyle="--")
    skewt.plot_lines_skewt(prof.pres, prof.tmpc, color="red", linewidth=2)
    skewt.plot_lines_skewt(
        prof.mlpcl.ptrace, prof.mlpcl.ttrace, color="darkred", linewidth=1,
        linestyle="--")
    skewt.plot_surface_values(prof.pres[0], prof.tmpc[0], prof.dwpc[0])
    skewt.plot_heights(prof.pres, prof.hght)
    skewt.plot_inflow(
        prof.ebottom, prof.etop, prof.ebotm, prof.etopm, prof.esrh[0])
    skewt.plot_max_lr(
        prof.max_lapse_rate_2_6[1], prof.max_lapse_rate_2_6[2],
        prof.max_lapse_rate_2_6[0])
    skewt.plot_sounding_levels(
        prof.mlpcl.lclpres, prof.mlpcl.lfcpres, prof.mlpcl.elpres)

    # barbs
    skewt.plot_barbs(prof.pres, prof.u, prof.v)

    # hodo
    skewt.plot_hodo(prof.hght / 1000, prof.u, prof.v)
    skewt.plot_storm_movers(
        prof.bunkers[0], prof.bunkers[1], prof.bunkers[2], prof.bunkers[3],
        *utils.vec2comp(prof.mean_lcl_el[0], prof.mean_lcl_el[1]))
    ebot_wind = interp.components(prof, prof.ebottom)
    etop_wind = interp.components(prof, prof.etop)
    skewt.plot_hodo_inflow(
        ebot_wind[0], ebot_wind[1], etop_wind[0], etop_wind[1], prof.bunkers[0],
        prof.bunkers[1])
    skewt.plot_critical_angle(prof.critical_angle)

    # theta
    skewt.plot_thetae(prof.pres, prof.thetae)

    # indicies
    sfc_1km_shear = utils.mag(prof.sfc_1km_shear[0], prof.sfc_1km_shear[1])
    sfc_3km_shear = utils.mag(prof.sfc_3km_shear[0], prof.sfc_3km_shear[1])
    sfc_6km_shear = utils.mag(prof.sfc_6km_shear[0], prof.sfc_6km_shear[1])
    eff_shear = utils.mag(prof.eff_shear[0], prof.eff_shear[1])
    mean_eff = utils.comp2vec(prof.mean_eff[0], prof.mean_eff[1])
    srw_eff = utils.comp2vec(prof.srw_eff[0], prof.srw_eff[1])
    bunkers_right = utils.comp2vec(prof.bunkers[0], prof.bunkers[1])
    bunkers_left = utils.comp2vec(prof.bunkers[2], prof.bunkers[3])
    wbz = interp.hght(prof, params.temp_lvl(prof, 0, wetbulb=True))
    fzl = interp.hght(prof, params.temp_lvl(prof, 0))
    indices = [
        ["PARCEL", "CAPE",                            "CIN",                              "LCL (m)",                           "LI",                            "LFC (m)",                           "EL (m)"],
        ["SFC",    utils.INT2STR(prof.sfcpcl.bplus),  utils.INT2STR(prof.sfcpcl.bminus),  utils.INT2STR(prof.sfcpcl.lclhght),  utils.INT2STR(prof.sfcpcl.li5),  utils.INT2STR(prof.sfcpcl.lfchght),  utils.INT2STR(prof.sfcpcl.elhght)],
        ["ML",     utils.INT2STR(prof.mlpcl.bplus),   utils.INT2STR(prof.mlpcl.bminus),   utils.INT2STR(prof.mlpcl.lclhght),   utils.INT2STR(prof.mlpcl.li5),   utils.INT2STR(prof.mlpcl.lfchght),   utils.INT2STR(prof.mlpcl.elhght)],
        ["FCST",   utils.INT2STR(prof.fcstpcl.bplus), utils.INT2STR(prof.fcstpcl.bminus), utils.INT2STR(prof.fcstpcl.lclhght), utils.INT2STR(prof.fcstpcl.li5), utils.INT2STR(prof.fcstpcl.lfchght), utils.INT2STR(prof.fcstpcl.elhght)],
        ["MU",     utils.INT2STR(prof.mupcl.bplus),   utils.INT2STR(prof.mupcl.bminus),   utils.INT2STR(prof.mupcl.lclhght),   utils.INT2STR(prof.mupcl.li5),   utils.INT2STR(prof.mupcl.lfchght),   utils.INT2STR(prof.mupcl.elhght)],
        ["", "", "", "", "", "", ""],
        ["",           "SRH",                         "Shear(kt)   ",               "MnWind(kt)",                                                           "  SRW(kt)", "", ""],
        ["SFC-1km",    utils.INT2STR(prof.srh1km[0]), utils.INT2STR(sfc_1km_shear), f"{utils.INT2STR(prof.mean_1km[0])}/{utils.INT2STR(prof.mean_1km[1])}", f"{utils.INT2STR(prof.srw_1km[0])}/{utils.INT2STR(prof.srw_1km[1])}", "", ""],
        ["SFC-3km",    utils.INT2STR(prof.srh3km[0]), utils.INT2STR(sfc_3km_shear), f"{utils.INT2STR(prof.mean_3km[0])}/{utils.INT2STR(prof.mean_3km[1])}", f"{utils.INT2STR(prof.srw_3km[0])}/{utils.INT2STR(prof.srw_3km[1])}", "", ""],
        ["Eff Inflow", utils.INT2STR(prof.esrh[0]),   utils.INT2STR(eff_shear),     f"{utils.INT2STR(mean_eff[0])}/{utils.INT2STR(mean_eff[1])}",           f"{utils.INT2STR(srw_eff[0])}/{utils.INT2STR(srw_eff[1])}", "", ""],
        ["SFC-6km",    "",                            utils.INT2STR(sfc_6km_shear), f"{utils.INT2STR(prof.mean_6km[0])}/{utils.INT2STR(prof.mean_6km[1])}", f"{utils.INT2STR(prof.srw_6km[0])}/{utils.INT2STR(prof.srw_6km[1])}", "1km red", "6km blue"],
        ["", "", "", "", "", "", ""],
        ["PWAT",   f"{utils.FLOAT2STR(prof.pwat, 2)}in", "3CAPE", utils.INT2STR(prof.mlpcl.b3km),                            "ConvT",     f"{utils.INT2STR(prof.convT)}°F", ""],
        ["K",      utils.INT2STR(prof.k_idx),            "DCAPE", utils.INT2STR(prof.dcape),                                 "SCP",       utils.FLOAT2STR(prof.scp, 1), ""],
        ["Mid RH", f"{utils.INT2STR(prof.mid_rh)}%",     "DownT", f"{utils.INT2STR(thermo.ctof(prof.dpcl_ttrace.max()))}°F", "STP (cin)", utils.FLOAT2STR(prof.stp_cin, 1), ""],
        ["Low RH", f"{utils.INT2STR(prof.low_rh)}%",     "MeanW", f"{utils.FLOAT2STR(prof.mean_mixr, 1)}g kg$^{{{-1}}}$",    "SHIP",      utils.FLOAT2STR(prof.ship, 1), ""],
        ["", "", "", "", "", "", ""],
        ["SFC-3km Γ", utils.FLOAT2STR(prof.lapserate_3km, 1),     "Bunkers R",  f"{utils.INT2STR(bunkers_right[0])}/{utils.INT2STR(bunkers_right[1])}kt", "WBZ",  f"{utils.INT2STR(wbz)}m", ""],
        ["3-6km Γ",   utils.FLOAT2STR(prof.lapserate_3_6km, 1),   "Bunkers L",  f"{utils.INT2STR(bunkers_left[0])}/{utils.INT2STR(bunkers_left[1])}kt",   "FZL",  f"{utils.INT2STR(fzl)}m", ""],
        ["850-500 Γ", utils.FLOAT2STR(prof.lapserate_850_500, 1), "Hazzard",     prof.watch_type,                                                         "MaxT", f"{utils.INT2STR(prof.maxT)}°F", ""],
        ["700-500 Γ", utils.FLOAT2STR(prof.lapserate_700_500, 1), "Precip Type", prof.precip_type,                                                        "",     "",                              ""],
    ]
    skewt.plot_indices(indices)

    wind_1km = interp.components(prof, interp.pres(prof, 1000))
    wind_6km = interp.components(prof, interp.pres(prof, 6000))
    skewt.plot_1_6km_barbs(wind_1km[0], wind_1km[1], wind_6km[0], wind_6km[1])

    skewt.save_image(
        f"../images_wrf/{init_time:%Y%m%d%H}/{ens}/skewt/{station_name}/",
        f"skewt_{station_name}_{fcst_time:%Y%m%d%H%M}")


def plot_hour(init_time, fcst_time, domain, ens):
    print("Reading data for", fcst_time, "domain", domain, flush=True)
    ds_wrf = read_data_wrf_sounding(init_time, fcst_time, domain, ens)
    ds_wrf.ds = add_cf_wrf(ds_wrf.ds)
    for station in wrf_sounding_stations[domain].values():
        plot_sounding(
            data=ds_wrf.ds,
            lon=station.lon,
            lat=station.lat,
            station_name=station.name,
            fcst_time=fcst_time,
            init_time=start_date,
            ens=ens)


def main():
    with mp.Pool(processes=24, maxtasksperchild=1) as pool:
        for domain in np.arange(
                config.wrf_max_domain, config.wrf_min_domain - 1, -1):
            for time in pd.date_range(
                    start_date, end_date, freq=config.wrf_freq[domain]):
                print(time, flush=True)
                pool.apply_async(
                    plot_hour,
                    args=(start_date, time, domain, args.ens, ),
                    error_callback=error_handler)
        pool.close()
        pool.join()


if __name__ == "__main__":
    mp.set_start_method("fork")
    main()
