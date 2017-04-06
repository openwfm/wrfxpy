#
# Dalton Burke
#


# Test correct functionality of the retrieval of level0 files

import subprocess
import datetime
import shutil

local_path = 'l0_test_ingest'

# Remove data from old tests
shutil.rmtree(local_path, ignore_errors=True)

current_time = datetime.datetime.utcnow()

ten_hours_ago = str(current_time - datetime.timedelta(hours=10)).replace(' ', '_')
five_hours_ago = str(current_time - datetime.timedelta(hours=5)).replace(' ', '_')
current_time = str(current_time).replace(' ', '_')

source_types = ['MODIS_AQUA', 'MODIS_TERRA', 'VIIRS_NPP']

# -----------------------------------------------------------------------
# Download all data sources from the last 5 hours

print "TESTING SOURCES FOR FILES IN LAST 5 HOURS\n"

for t in source_types:
    print "\nRETRIEVING %s FILES FROM THE LAST 5 HOURS WITH CALL:" % t
    print './level0_retr.sh %s %s %s %s \n' % (t, five_hours_ago, current_time, local_path)

    subprocess.call(['./level0_retr.sh', t, five_hours_ago, current_time, local_path])

print "\nDONE RETRIEVING FILES FROM LAST 5 HOURS \n\n"

# -----------------------------------------------------------------------
# Download all data sources from the last 10 hours
# (some data we should already have, so those should be skipped)

print "TESTING SOURCES FOR FILES IN LAST 10 HOURS\n"

for t in source_types:
    print "\nRETRIEVING %s FILES FROM THE LAST 10 HOURS WITH CALL:" % t
    print './level0_retr.sh %s %s %s %s \n' % (t, ten_hours_ago, current_time, local_path)

    subprocess.call(['./level0_retr.sh', t, ten_hours_ago, current_time, local_path])

print "\nDONE RETRIEVING FILES FROM LAST 10 HOURS"
