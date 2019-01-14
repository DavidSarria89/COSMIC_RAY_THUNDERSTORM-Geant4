#
# Make example source for MCNP and MCNPX
#
# convert multiple blank lines to one using sed or awk
# http://www.zazzybob.com/bin/blankcat.sed.html

default: source.F ../mcnp/source.F90

# MCNP
../mcnp/source.F90: source.master
	cpp -P $< $@

# MCNPX
source.F: source.master
	cpp -P -D=MCNPX $< $@
