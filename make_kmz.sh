 #!/usr/bin/env bash
cd $(dirname "$0")
./ssh_command.sh wrfxweb/make_kmz.sh $*
