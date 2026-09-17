# -*- makefile-gmake -*-
# Generated automatically! Do not edit

.PHONY:
all: pkpm pkpm-unit pkpm-regression ## Build only specified Apps (pkpm)

.PHONY: install
install: pkpm-install ## Install gkeyll executable
	cd gkeyll && ${MAKE} -f Makefile-gkeyll install

