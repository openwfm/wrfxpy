#
# Dalton Burke
#


# Test correct functionality of the retrieval of level0 files

import subprocess
import shlex
import datetime


current_time = datetime.datetime.utcnow()
ten_hours_ago = str(current_time - datetime.timedelta(hours=10)).replace(' ', '_')
five_hours_ago = str(current_time - datetime.timedelta(hours=5)).replace(' ', '_')
current_time = str(current_time).replace(' ', '_')

print "\nRETRIEVING AQUA FILES FROM THE LAST 5 HOURS WITH CALL:"
print ('./level0_retr.sh %s %s %s %s \n' %
       ('MODIS_AQUA', five_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'MODIS_AQUA', five_hours_ago, current_time, '~/level0test/ingest'])

print "\nRETRIEVING TERRA FILES FROM THE LAST 5 HOURS WITH CALL"
print ('./level0_retr.sh %s %s %s %s\n' %
       ('MODIS_TERRA', five_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'MODIS_TERRA', five_hours_ago, current_time, '~/level0test/ingest'])

print "\nRETRIEVING VIIRS FILES FROM THE LAST 5 HOURS WITH CALL:"
print ('./level0_retr.sh %s %s %s %s' %
       ('VIIRS_NPP', five_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'VIIRS_NPP', five_hours_ago, current_time, '~/level0test/ingest'])

print "\nDONE RETRIEVING FILES FROM LAST 5 HOURS \n\n"

print "\nRETRIEVING AQUA FILES FROM LAST 10 HOURS (We should already have some) WITH CALL:"
print ('./level0_retr.sh %s %s %s %s\n' %
       ('MODIS_AQUA', ten_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'MODIS_AQUA', ten_hours_ago, current_time, '~/level0test/ingest'])

print "\nRETRIEVING TERRA FILES FROM LAST 10 HOURS (We should already have some) WITH CALL:"
print ('./level0_retr.sh %s %s %s %s\n' %
       ('MODIS_TERRA', ten_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'MODIS_TERRA', ten_hours_ago, current_time, '~/level0test/ingest'])

print "\nRETRIEVING VIIRS FILES FROM LAST 10 HOURS (We should already have some) WITH CALL:"
print ('./level0_retr.sh %s %s %s %s\n' %
       ('VIIRS_NPP', ten_hours_ago, current_time, '~/level0test/ingest'))

subprocess.call(['./level0_retr.sh', 'VIIRS_NPP', ten_hours_ago, current_time, '~/level0test/ingest'])

print "\nDONE RETRIEVING FILES FROM LAST 10 HOURS"
