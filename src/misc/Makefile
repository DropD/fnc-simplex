#######################################
# Redirective Makefile, for free use
# In case you want to use a centralised
# Makefile in some other directory.
#
# Use
#  DEFAULT_TARGET=
# to use the target Makefile default


TARGETDIR=".."
CUSTOM_MAKEFLAGS="-j6"
DEFAULT_TARGET=


########################################################################

.PHONY: all what

all:
	@cd $(TARGETDIR) && \
	`which time` -f "\treal: %e \tuser: %U \tsys: %S  \twall clock: %e" \
	make $(CUSTOM_MAKEFLAGS) $(DEFAULT_TARGET)

%:
	@cd $(TARGETDIR) && \
	`which time` -f "\treal: %e \tuser: %U \tsys: %S  \twall clock: %e" \
	make $(CUSTOM_MAKEFLAGS) $@

what:
	@echo "AVAILABLE TARGETS:"
	@cat $(TARGETDIR)"/Makefile" | grep -Pv '^(#|\s)' \
		| grep -v '^.*:=' \
		| grep --only-matching '^.*:' \
		| sed 's/://' | grep -v ',' \
		| grep -v '^.PHONY' | grep -v '^.SUFFIXES'
