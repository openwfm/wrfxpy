 #!/usr/bin/env bash

#get CFSR pressure data
echo "getting CFSR pressure data"
cd /share_home/akochans/wrfxpy_new
./grib_retr.sh CFSR_P 2015-09-25_00:00:00 2015-09-29_00:00:00 ingest

cd  /share_home/akochans/wrfxpy_new/wksp/wfc-Fishlake_5d_FSHU1-2016-09-11_00:00:00-48/wps

rm GRIBFILE.*
rm COLMET_*

#link namelist wps with ungrib prefix set to COLMET_P
rm namelist.wps
ln -s namelist.wps_P namelist.wps

#link Variable Table for CFSR pressure data
rm Vtable
ln -s /share_home/akochans/wrfxpy_new/etc/vtables/Vtable.CFSR_press_pgbh06 ./Vtable 

#link CFSR_P grib files
echo "linking CFSR pressure grib files"
./link_grib.csh ../../../ingest/2016/201609/20160910/*pgrbh00.grib2 ../../../ingest/2016/201609/20160911/*pgrbh00.grib2 ../../../ingest/2016/201609/20160912/*pgrbh00.grib2 ../../../ingest/2016/201609/20160913/*pgrbh00.grib2


#ungrib files to generate COLMET_P intermediate files
echo "ungribing pressure files"
./ungrib.exe

#at that point the CFSR pressure data should be processed to the wps intermediate format, and a set of COLMET_P:XXXXX files should be generated
#the next step is to do the same, but with the CFSR surface data
#---------------------------------------------------------------

#get CFSR surface data
echo "getting CFSR surface data"
cd /share_home/akochans/wrfxpy_new
./grib_retr.sh CFSR_S 2016-09-10_00:00:00 2016-09-12_00:00:00 ingest

cd  /share_home/akochans/wrfxpy_new/wksp/wfc-Fishlake_5d_FSHU1-2016-09-11_00:00:00-48/wps

#link namelist wps with ungrib prefix set to COLMET_S (not to overwrite the previously generated COLMET_P pressure files)
rm namelist.wps
ln -s namelist.wps_S namelist.wps

#link Variable Table for CFSR surface data
rm Vtable
ln -s /share_home/akochans/wrfxpy_new/etc/vtables/Vtable.CFSR_sfc_flxf06 ./Vtable

echo "linking CFSR surface data"
#link CFSR_S grib files
./link_grib.csh ../../../ingest/2016/201609/20160910/*sfluxgrbf00.grib2 ../../../ingest/2016/201609/20160911/*sfluxgrbf00.grib2 ../../../ingest/2016/201609/20160912/*sfluxgrbf00.grib2 ../../../ingest/2016/201609/20160913/*sfluxgrbf00.grib2

#ungrib files to generate COLMET_S intermediate files
echo "ungribing surface files"
./ungrib.exe

#at that point CFSR surface data should be processed to the wps intermediate format, and a set of COLMET_S:XXXXX files should be generated
#we have both pressure and surface data processed (COLET_P, and COLMET_S files), so we are ready to run metgrid and then real as usual.
#--------------------------------------------------------------------------------------------------------------------------------------



echo "running metgrid.exe"
#heaving COLMET_P and COLMET_S internediate files run metgrid (there are two prefixes in the mtgrid section of namelist.wps, so both sets of files get processed)
./metgrid.exe

# after metgrid is done link met_em files to wrf directory and run real.exe, anc then wrf.exe

echo "running real.exe"
cd ../wrf/
ln -s ../wps/met_em* ./

./real.exe 

echo "done! ready to start wrf.exe"
