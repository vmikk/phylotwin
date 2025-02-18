#!/bin/bash

## Quote species names in a Nexus file

sed -E "/TAXLABELS/,/;/ s/^\t\t(.+)$/\t\t'\1'/; /TRANSLATE/,/;/ s/^\t\t([0-9]+\t)([^,]+)(,?)$/\t\t\1'\2'\3/" "$1"

