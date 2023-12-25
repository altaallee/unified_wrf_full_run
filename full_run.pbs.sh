#!/bin/bash

uw=true
hrrr=true
conus=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--date)
            date="$2"
            shift
            shift
            ;;
        -r|--run)
            run="$2"
            shift
            shift
            ;;
        -t|--total_hours)
            total_hours="$2"
            shift
            shift
            ;;
        -f|--forward_hours)
            forward_hours="$2"
            shift
            shift
            ;;
        -*|--*)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

if [ -z "${date}" ]; then
    echo Specify GFS init date with --date YYYYMMDD
    exit 1
fi
if [ -z "${run}" ]; then
    echo Specify GFS init run with --run HH
    exit 1
fi
if [ -z "${total_hours}" ]; then
    total_hours=72
    echo "Total WRF run hours set to ${total_hours}"
    echo Specify total_hours with --total_hours
fi
if [ -z "${forward_hours}" ]; then
    forward_hours=0
    echo "Forward hours set to ${forward_hours}"
    echo Specify forward_hours with --forward_hours
fi

wrf_date=$(date -d "$date $run+${forward_hours}hours" "+%Y%m%d")
wrf_run=$(date -d "$date $run+${forward_hours}hours" "+%H")

echo "GFS ICBC date" $date $run
echo "WRF forecast start date" $wrf_date $wrf_run

rm files_page.txt
rm pbs_out_metgrid
rm pbs_out_wrf_conus
rm pbs_out_wrf_hrrr
rm pbs_out_wrf_uw
rm pbs_out_wrf_postprocess_conus
rm pbs_out_wrf_postprocess_hrrr
rm pbs_out_wrf_postprocess_uw
rm -rf ../gfs_ICBC_wrf
rm ../WPS/GRIBFILE*
rm ../WPS/FILE:*
rm ../WPS/geogrid.log.*
rm ../WPS/met_em.*
rm ../WPS/metgrid.log.*
rm ../WPS/out_geogrid
rm ../WPS/out_metgrid
rm ../WPS/out_ungrib

rm ../WRF/run_uw/met_em.*
rm ../WRF/run_uw/out_real
rm ../WRF/run_uw/out_wrf
rm ../WRF/run_uw/rsl.error.*
rm ../WRF/run_uw/rsl.out.*
rm ../WRF/run_uw/wrfbdy*
rm ../WRF/run_uw/wrfinput*
rm ../WRF/run_uw/wrfrst*

rm ../WRF/run_hrrr/met_em.*
rm ../WRF/run_hrrr/out_real
rm ../WRF/run_hrrr/out_wrf
rm ../WRF/run_hrrr/rsl.error.*
rm ../WRF/run_hrrr/rsl.out.*
rm ../WRF/run_hrrr/wrfbdy*
rm ../WRF/run_hrrr/wrfinput*
rm ../WRF/run_hrrr/wrfrst*

rm ../WRF/run_conus/met_em.*
rm ../WRF/run_conus/out_real
rm ../WRF/run_conus/out_wrf
rm ../WRF/run_conus/rsl.error.*
rm ../WRF/run_conus/rsl.out.*
rm ../WRF/run_conus/wrfbdy*
rm ../WRF/run_conus/wrfinput*
rm ../WRF/run_conus/wrfrst*

mkdir ../gfs_ICBC_wrf
curl "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${date}/${run}/atmos/" > files_page.txt
for hour in $( eval echo {0..$total_hours..1} ); do
    pad_hour=$(printf "%03d\\n" $(( $hour + $forward_hours )))
    echo $pad_hour
    for recheck in {30..0}; do
        if grep "gfs.t${run}z.pgrb2.0p25.f${pad_hour}" files_page.txt; then
            curl "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t${run}z.pgrb2.0p25.f${pad_hour}&lev_0-0.1_m_below_ground=on&lev_0.1-0.4_m_below_ground=on&lev_0.4-1_m_below_ground=on&lev_1000_mb=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_10_mb=on&lev_1-2_m_below_ground=on&lev_150_mb=on&lev_15_mb=on&lev_200_mb=on&lev_20_mb=on&lev_250_mb=on&lev_2_m_above_ground=on&lev_300_mb=on&lev_30_mb=on&lev_350_mb=on&lev_3658_m_above_mean_sea_level=on&lev_400_mb=on&lev_40_mb=on&lev_450_mb=on&lev_500_mb=on&lev_50_mb=on&lev_550_mb=on&lev_600_mb=on&lev_650_mb=on&lev_700_mb=on&lev_70_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_975_mb=on&lev_max_wind=on&lev_mean_sea_level=on&lev_surface=on&lev_tropopause=on&var_HGT=on&var_LAND=on&var_MSLET=on&var_PRES=on&var_PRMSL=on&var_RH=on&var_SNOD=on&var_SOILW=on&var_TMP=on&var_TSOIL=on&var_UGRD=on&var_VGRD=on&var_WEASD=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.${date}%2F${run}%2Fatmos" > "../gfs_ICBC_wrf/gfs.${date}.${run}z.f${pad_hour}"
            break
        fi
		if [ $recheck -ne 0 ]; then
            echo "GFS file not found."
            sleep 60
            curl "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${date}/${run}/atmos/" > files_page.txt
			continue
		else
			exit
		fi
    done
    sleep 1
done

if [ $uw = true ]; then
    python3 update_namelist.py --date $wrf_date$wrf_run --hours $total_hours --ens uw
    mkdir -p ../wrfout/$wrf_date$wrf_run/uw
fi
if [ $hrrr = true ]; then
    python3 update_namelist.py --date $wrf_date$wrf_run --hours $total_hours --ens hrrr
    mkdir -p ../wrfout/$wrf_date$wrf_run/hrrr
fi
if [ $conus = true ]; then
    python3 update_namelist.py --date $wrf_date$wrf_run --hours $total_hours --ens conus
    mkdir -p ../wrfout/$wrf_date$wrf_run/conus
fi

cd ../WPS
./link_grib.csh ../gfs_ICBC_wrf/gfs.${date}.${run}z*
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable
./ungrib.exe > out_ungrib

cd ../wrf_full_run

metgrid_id=$(qsub run_metgrid.pbs.sh)

if [ $uw = true ]; then
    wrf_uw_id=$(qsub -W depend=afterok:$metgrid_id run_wrf_uw.pbs.sh)
    maps_uw_id=$(qsub -W depend=afterok:$wrf_uw_id -N wrf_postprocess_1 -o pbs_out_wrf_postprocess_uw -S /bin/bash -l select=1:ncpus=24:model=has -l walltime=1:05:00 -j oe -W group_list=s2395 -m abe -v DATE=$wrf_date,RUN=$wrf_run,HOURS=$total_hours,ENS=uw $(pwd)/run_postprocess.pbs.sh)
fi

if [ $hrrr = true ]; then
    wrf_hrrr_id=$(qsub -W depend=afterok:$metgrid_id run_wrf_hrrr.pbs.sh)
    maps_hrrr_id=$(qsub -W depend=afterok:$wrf_hrrr_id -N wrf_postprocess_2 -o pbs_out_wrf_postprocess_hrrr -S /bin/bash -l select=1:ncpus=24:model=has -l walltime=1:05:00 -j oe -W group_list=s2395 -m abe -v DATE=$wrf_date,RUN=$wrf_run,HOURS=$total_hours,ENS=hrrr $(pwd)/run_postprocess.pbs.sh)
fi

if [ $conus = true ]; then
    wrf_conus_id=$(qsub -W depend=afterok:$metgrid_id run_wrf_conus.pbs.sh)
    maps_conus_id=$(qsub -W depend=afterok:$wrf_conus_id -N wrf_postprocess_3 -o pbs_out_wrf_postprocess_conus -S /bin/bash -l select=1:ncpus=24:model=has -l walltime=1:05:00 -j oe -W group_list=s2395 -m abe -v DATE=$wrf_date,RUN=$wrf_run,HOURS=$total_hours,ENS=conus $(pwd)/run_postprocess.pbs.sh)
fi

qstat -u alee31
