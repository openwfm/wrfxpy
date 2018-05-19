import logging
from utils import esmf_to_utc, utc_to_esmf, round_time_to_hour, timedelta_hours, load_sys_cfg
from ingest.grib_source import NAM218

# ESMF date: YYYY-MM-DD_hh:mm:ss

def test_time():

    cfg = load_sys_cfg()

    g = NAM218(cfg)

    cycle_start_esmf = "2005-01-05_00:00:00"
    cycle_start =  esmf_to_utc(cycle_start_esmf)
    from_utc = esmf_to_utc("2005-01-05_00:00:00")
    to_utc = esmf_to_utc("2005-01-06_00:00:01")



    print cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour

    fc_hours = int((to_utc - cycle_start).total_seconds())/3600


    delta = from_utc - cycle_start
    print str(delta)
    print delta.days, delta.seconds, delta.total_seconds()
    print (to_utc - cycle_start).total_seconds() 
    print timedelta_hours(to_utc - cycle_start)
    print timedelta_hours(to_utc - cycle_start, False)

    g.forecast_times(cycle_start, from_utc, to_utc)

     
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


test_time()
