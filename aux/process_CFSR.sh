 #!/usr/bin/env bash
cd $(dirname "$0")/..
wrfxpy=`pwd`
echo wrfxpy=$wrfxpy
tmp=$wrfxpy/wksp/tmp
echo tmp=$tmp
mkdir -p $tmp

wps=`awk '/wps_install/ {print $2}' $wrfxpy/etc/conf.json | sed 's/"//g' | sed 's/,//'`
echo wps=$wps

cd $tmp
rm -f GRIBFILE.*
rm -f COLMET_*

echo get CFSR pressure data
$wrfxpy/grib_retr.sh CFSR_P 2016-09-10_00:00:00 2016-09-13_00:00:00 ingest

echo link namelist wps with ungrib prefix set to COLMET_P
rm -f namelist.wps
ln -s $wrfxpy/etc/nlists/default.wps_P namelist.wps
ls -l namelist.wps
#cat namelist.wps

echo link Variable Table for CFSR pressure data
rm -f Vtable
ln -s $wrfxpy/etc/vtables/Vtable.CFSR_press_pgbh06 ./Vtable 
ls -l Vtable
#cat Vtable

#link CFSR_P grib files
echo linking CFSR pressure grib files
ls -l $wrfxpy/ingest/2016/201609/20160910/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160911/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160912/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160913/*pgrbh00.grib2
$wps/link_grib.csh $wrfxpy/ingest/2016/201609/20160910/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160911/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160912/*pgrbh00.grib2 $wrfxpy/ingest/2016/201609/20160913/*pgrbh00.grib2

#ungrib files to generate COLMET_P intermediate files
echo "ungribing pressure files"
$wps/ungrib.exe

#at that point the CFSR pressure data should be processed to the wps intermediate format, and a set of COLMET_P:XXXXX files should be generated
#the next step is to do the same, but with the CFSR surface data
#---------------------------------------------------------------

#get CFSR surface data
echo "getting CFSR surface data"
$wrfxpy/grib_retr.sh CFSR_S 2016-09-10_00:00:00 2016-09-13_00:00:00 ingest

#link namelist wps with ungrib prefix set to COLMET_S (not to overwrite the previously generated COLMET_P pressure files)
rm -f namelist.wps
ln -s $wrfxpy/etc/nlists/default.wps_S namelist.wps
ls -l namelist.wps

#link Variable Table for CFSR surface data
rm Vtable
ln -s $wrfxpy/etc/vtables/Vtable.CFSR_sfc_flxf06 ./Vtable

echo "linking CFSR surface data"
#link CFSR_S grib files
$wps/link_grib.csh $wrfxpy/ingest/2016/201609/20160910/*sfluxgrbf00.grib2 $wrfxpy/ingest/2016/201609/20160911/*sfluxgrbf00.grib2 $wrfxpy/ingest/2016/201609/20160912/*sfluxgrbf00.grib2 $wrfxpy/ingest/2016/201609/20160913/*sfluxgrbf00.grib2

#ungrib files to generate COLMET_S intermediate files
echo "ungribing surface files"
$wps/ungrib.exe

#at that point CFSR surface data should be processed to the wps intermediate format, and a set of COLMET_S:XXXXX files should be generated
#we have both pressure and surface data processed (COLET_P, and COLMET_S files), so we are ready to run metgrid and then real as usual.
#--------------------------------------------------------------------------------------------------------------------------------------



echo "running metgrid.exe"
#heaving COLMET_P and COLMET_S internediate files run metgrid (there are two prefixes in the mtgrid section of namelist.wps, so both sets of files get processed)
$wps/metgrid.exe

# after metgrid is done link met_em files to wrf directory and run real.exe, anc then wrf.exe

#echo "running real.exe"

#./real.exe 

#echo "done! ready to start wrf.exe"
