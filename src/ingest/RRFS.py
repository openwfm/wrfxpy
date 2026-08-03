from ingest.grib_forecast import GribForecast


class RRFS(GribForecast):
    """
    RRFSv1 forecast grib source from NOAA AWS.

    RRFS is split into pressure-level and 2D/surface products. This container
    class holds the common metadata and WPS configuration, while subclasses
    specify the domain and product type.

    For WPS, use the pressure product together with the 2D/surface product:

        prslev: 3D pressure-level atmospheric fields
        2dfld : surface, near-surface, land, soil, snow, radiation, precipitation fields
    """

    domain = None
    grid = None
    product = None
    product_root = "rrfs_public"

    def __init__(self, arg):
        super(RRFS, self).__init__(arg)

    def vtables(self):
        """
        Returns the variable tables that must be linked for RRFS.

        There is no standard WPS RRFS Vtable in many WPS installations.
        The HRRR Vtable/METGRID pair is the best starting point because RRFS
        prslev + 2dfld follows a similar WRF-oriented split between pressure
        fields and 2D/surface fields.

        This should be validated with:
            strings COLMET_P:YYYY-MM-DD_HH | grep -Ei "TT|UU|VV|RH|HGT"
            strings COLMET_S:YYYY-MM-DD_HH | grep -Ei "SOILT|SOILM|PSFC|SKINTEMP"
            ncdump -h met_em.d01.YYYY-MM-DD_HH:00:00.nc | grep -E "NUM_METGRID|BOTTOM-TOP"
        """
        return {
            'geogrid_vtable': 'GEOGRID.TBL',
            'ungrib_vtable': 'Vtable.RRFS',
            'metgrid_vtable': 'METGRID.TBL.RRFS',
        }

    def namelist_keys(self):
        """
        Returns namelist.input keys for RRFS.

        The values below assume HRRR-like pressure-level output and 9 soil levels.
        Confirm from met_em files after the first successful metgrid run.
        """
        return {
            'domains': {
                'num_metgrid_levels': 46,
                'num_metgrid_soil_levels': 9,
            }
        }

    def namelist_wps_keys(self):
        """
        Base container does not define a WPS prefix. Product subclasses do.
        """
        return None

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required RRFS GRIB2 files.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of forecast hours
        """
        if self.domain is None or self.grid is None or self.product is None:
            raise NotImplementedError(
                "RRFS subclasses must define domain, grid, and product."
            )

        path_tmpl = (
            '{root}/rrfs.%04d%02d%02d/%02d/'
            'rrfs.t%02dz.{product}.{grid}.f%03d.{domain}.grib2'
        )

        grib_files = [
            path_tmpl.format(
                root=self.product_root,
                product=self.product,
                grid=self.grid,
                domain=self.domain,
            ) % (
                cycle_start.year,
                cycle_start.month,
                cycle_start.day,
                cycle_start.hour,
                cycle_start.hour,
                fh,
            )
            for fh in fc_list
        ]

        return grib_files

    id = "RRFS"
    info_url = "https://rapidrefresh.noaa.gov/RRFS/"
    info_aws = "https://registry.opendata.aws/noaa-rrfs/"
    info_text = "NOAA RRFSv1 Rapid Refresh Forecast System"
    info = "The Rapid Refresh Forecast System (RRFSv1)"
    remote_url = ["s3://noaa-rrfs-pds/"]
    browse_aws = "https://noaa-rrfs-pds.s3.amazonaws.com/index.html"


class RRFS_CONUS(RRFS):
    """
    RRFSv1 CONUS container.

    The CONUS products use the 3-km CONUS grid:

        rrfs.tCCz.prslev.3km.fFFF.conus.grib2
        rrfs.tCCz.2dfld.3km.fFFF.conus.grib2
    """

    def __init__(self, arg):
        super(RRFS_CONUS, self).__init__(arg)

    id = "RRFS_CONUS"
    domain = "conus"
    grid = "3km"
    info_text = "NOAA RRFSv1 CONUS 3-km Forecast"
    cycle_hours = 3
    period_hours = 1
    hours_behind_real_time = 1
    grib_forecast_hours_periods = [{'hours': 84, 'period': 1}]


class RRFS_CONUS_P(RRFS_CONUS):
    """
    RRFSv1 CONUS pressure-level product.

    Provides 3D pressure-level atmospheric fields.
    """

    def __init__(self, arg):
        super(RRFS_CONUS_P, self).__init__(arg)

    def namelist_wps_keys(self):
        return {
            'ungrib': {'prefix': 'COLMET_P'},
            'metgrid': {'fg_name': ['COLMET_S', 'COLMET_P']},
        }

    id = "RRFS_CONUS_P"
    product = "prslev"
    prefix = "COLMET_P"


class RRFS_CONUS_S(RRFS_CONUS):
    """
    RRFSv1 CONUS 2D/surface product.

    Provides surface, near-surface, land, soil, snow, radiation,
    precipitation, and other 2D fields.
    """

    def __init__(self, arg):
        super(RRFS_CONUS_S, self).__init__(arg)

    def namelist_wps_keys(self):
        return {
            'ungrib': {'prefix': 'COLMET_S'},
            'metgrid': {'fg_name': ['COLMET_S', 'COLMET_P']},
        }

    id = "RRFS_CONUS_S"
    product = "2dfld"
    prefix = "COLMET_S"


class RRFS_NA(RRFS):
    """
    RRFSv1 North America container.
    """

    def __init__(self, arg):
        super(RRFS_NA, self).__init__(arg)

    id = "RRFS_NA"
    domain = "na"
    grid = "13km"
    info_text = "NOAA RRFSv1 North America Forecast"
    cycle_hours = 6
    period_hours = 3
    hours_behind_real_time = 1
    grib_forecast_hours_periods = [{'hours': 84, 'period': 3}]


class RRFS_NA_P(RRFS_NA):
    """
    RRFSv1 North America pressure-level product.
    """

    def __init__(self, arg):
        super(RRFS_NA_P, self).__init__(arg)

    def namelist_wps_keys(self):
        return {
            'ungrib': {'prefix': 'COLMET_P'},
            'metgrid': {'fg_name': ['COLMET_S', 'COLMET_P']},
        }

    id = "RRFS_NA_P"
    product = "prslev"
    prefix = "COLMET_P"


class RRFS_NA_S(RRFS_NA):
    """
    RRFSv1 North America 2D/surface product.
    """

    def __init__(self, arg):
        super(RRFS_NA_S, self).__init__(arg)

    def namelist_wps_keys(self):
        return {
            'ungrib': {'prefix': 'COLMET_S'},
            'metgrid': {'fg_name': ['COLMET_S', 'COLMET_P']},
        }

    id = "RRFS_NA_S"
    product = "2dfld"
    prefix = "COLMET_S"