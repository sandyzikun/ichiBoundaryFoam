# GNU Make front end for the OpenFOAM wmake build.

SHELL := /bin/sh

PROJECT := ichiBoundaryFoam
SOURCE_DIR := src
MARKER_TEST_DIR := tests/markerOperatorTest
UUID := B361474C_D507_11F1_9472_A3B0B9D8CE68
APP_BIN ?= $(CURDIR)/platforms/$(WM_OPTIONS)/bin
EXE_NAME ?= $(PROJECT)
EXE_PATH ?= $(APP_BIN)/$(EXE_NAME)
MARKER_TEST_EXE_NAME ?= ichiMarkerOperatorTest
MARKER_TEST_EXE_PATH ?= $(APP_BIN)/$(MARKER_TEST_EXE_NAME)
MARKER_SUBITER_EXE_NAME ?= ichiBoundaryFoam-markerSubIter
MARKER_SUBITER_EXE_PATH ?= $(APP_BIN)/$(MARKER_SUBITER_EXE_NAME)

## [INFO] Compiling the SDF-based direct-forcing IBM backend.
ENABLE_DIRECT_FORCING ?= 1

## [INFO] Compiling the unmodified interFoam execution path.
##  [NOTE] No geometry or IBM object is created in this mode.
PURE_INTERFOAM ?= 0

## [INFO] Compiling immersed-surface slip sampling and fields.
##  [NOTE] Disable this option to remove its run-time overhead.
ENABLE_SLIP_SAMPLING ?= 1

## [INFO] Compiling IBM constraint-impulse load accounting.
##  [NOTE] This is diagnostic only and does not feed loads back to
##  the prescribed or stationary immersed body.
ENABLE_LOAD_CALCULATION ?= 1

## [INFO] Compiling the SDF-based immersed turbulence-wall model.
##  [NOTE] The model remains disabled unless it is also enabled in
##  constant/ichiBoundaryGeometry.
ENABLE_TURBULENCE_WALL_TREATMENT ?= 0

## [INFO] Solving a global pressure-correction equation after every
##  MDF forcing pass.
##  [NOTE] The default MDF path repeats forcing on the predictor and
##  uses only the regular pressure projection.
ENABLE_MDF_PRESSURE_PROJECTION ?= 0
MDF_PP := $(ENABLE_MDF_PRESSURE_PROJECTION)

## [INFO] Compiling the experimental IBM-aware projection prototype.
##  [NOTE] This only makes the run-time-selectable backend available.
##  The regular pressure projection remains the default path.
ENABLE_IBM_AWARE_PROJECTION ?= 0
IBM_AWARE_PP := $(ENABLE_IBM_AWARE_PROJECTION)

## [INFO] Compiling passive one-way kinematic particle tracking.
##  [NOTE] Particles sample the carrier flow but do not feed momentum,
##  volume fraction, or any other source back to the Eulerian fields.
ENABLE_ONE_WAY_DPM ?= 0

PROJECT_DEFINES := \
	-DICHIBOUNDARYFOAM_PROJECT_UUID_$(UUID)=1 \
	-DICHIBOUNDARYFOAM_ENABLE_DIRECT_FORCING=$(ENABLE_DIRECT_FORCING) \
	-DICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING=$(ENABLE_SLIP_SAMPLING) \
	-DICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION=$(ENABLE_LOAD_CALCULATION) \
	-DICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT=$(ENABLE_TURBULENCE_WALL_TREATMENT) \
	-DICHIBOUNDARYFOAM_ENABLE_MDF_PRESSURE_PROJECTION=$(MDF_PP) \
	-DICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION=$(IBM_AWARE_PP) \
	-DICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM=$(ENABLE_ONE_WAY_DPM) \
	-DICHIBOUNDARYFOAM_PURE_INTERFOAM=$(PURE_INTERFOAM)

BUILD_FLAGS := $(strip \
	$(FOAM_EXTRA_CXXFLAGS) $(PROJECT_DEFINES))

.PHONY: all build ibm ibm-aware marker-subiteration pure-interfoam marker-operator-test \
	clean-marker-operator-test debug clean distclean show-config help

all: build

build: check-openfoam
	@compiler="$$(mpicxx --showme:command)"; \
	clean_path="$$(printf '%s' "$$PATH" | tr ':' '\n' | \
		grep -v '/lib/spack/env' | paste -sd: -)"; \
	PATH="$${compiler%/*}:$$clean_path" \
	FOAM_USER_APPBIN='$(APP_BIN)' \
	ICHIBOUNDARYFOAM_EXE='$(EXE_PATH)' \
	FOAM_EXTRA_CXXFLAGS='$(BUILD_FLAGS)' \
		wmake $(SOURCE_DIR)

## [INFO] Rebuilding the regular IBM development solver.
ibm: clean
	@$(MAKE) build PURE_INTERFOAM=0

## [INFO] Rebuilding an IBM-aware projection test executable.
##  [NOTE] Use EXE_NAME or EXE_PATH to keep it separate from the
##  regular development solver.
ibm-aware: clean
	@$(MAKE) build PURE_INTERFOAM=0 \
		ENABLE_IBM_AWARE_PROJECTION=1

## [INFO] Rebuilding the Constant--Pinelli marker sub-iteration
##  executable under a persistent, non-production binary name.
marker-subiteration:
	@$(MAKE) ibm-aware \
		EXE_NAME='$(MARKER_SUBITER_EXE_NAME)' \
		EXE_PATH='$(MARKER_SUBITER_EXE_PATH)'

## [INFO] Rebuilding the pure interFoam regression solver.
pure-interfoam: clean
	@$(MAKE) build PURE_INTERFOAM=1 ENABLE_DIRECT_FORCING=0 \
		ENABLE_ONE_WAY_DPM=0

## [INFO] Building the standalone marker-operator verification utility.
##  [NOTE] This executable reads a mesh and analytic fields but solves
##  no transport, momentum, or pressure equation.
marker-operator-test: check-openfoam
	@compiler="$$(mpicxx --showme:command)"; \
	clean_path="$$(printf '%s' "$$PATH" | tr ':' '\n' | \
		grep -v '/lib/spack/env' | paste -sd: -)"; \
	PATH="$${compiler%/*}:$$clean_path" \
	FOAM_USER_APPBIN='$(APP_BIN)' \
	ICHI_MARKER_TEST_EXE='$(MARKER_TEST_EXE_PATH)' \
	FOAM_EXTRA_CXXFLAGS='$(BUILD_FLAGS)' \
		wmake $(MARKER_TEST_DIR)

clean-marker-operator-test: check-openfoam
	@ICHI_MARKER_TEST_EXE='$(MARKER_TEST_EXE_PATH)' \
		wclean $(MARKER_TEST_DIR)

debug: check-openfoam
	@compiler="$$(mpicxx --showme:command)"; \
	clean_path="$$(printf '%s' "$$PATH" | tr ':' '\n' | \
		grep -v '/lib/spack/env' | paste -sd: -)"; \
	PATH="$${compiler%/*}:$$clean_path" \
	FOAM_USER_APPBIN='$(APP_BIN)' \
	ICHIBOUNDARYFOAM_EXE='$(EXE_PATH)' \
	FOAM_EXTRA_CXXFLAGS='$(BUILD_FLAGS)' \
		wmake -debug $(SOURCE_DIR)

clean: check-openfoam
	@FOAM_USER_APPBIN='$(APP_BIN)' \
	ICHIBOUNDARYFOAM_EXE='$(EXE_PATH)' \
		wclean $(SOURCE_DIR)

distclean: clean
	@rm -f "$(EXE_PATH)"

show-config:
	@printf '%-20s %s\n' PROJECT $(PROJECT)
	@printf '%-20s %s\n' UUID $(UUID)
	@printf '%-20s %s\n' APP_BIN $(APP_BIN)
	@printf '%-20s %s\n' EXE_NAME $(EXE_NAME)
	@printf '%-20s %s\n' EXE_PATH $(EXE_PATH)
	@printf '%-20s %s\n' MARKER_TEST_EXE \
		$(MARKER_TEST_EXE_PATH)
	@printf '%-20s %s\n' MARKER_SUBITER_EXE \
		$(MARKER_SUBITER_EXE_PATH)
	@printf '%-20s %s\n' ENABLE_DIRECT_FORCING \
		$(ENABLE_DIRECT_FORCING)
	@printf '%-20s %s\n' PURE_INTERFOAM $(PURE_INTERFOAM)
	@printf '%-20s %s\n' ENABLE_SLIP_SAMPLING \
		$(ENABLE_SLIP_SAMPLING)
	@printf '%-20s %s\n' LOAD_CALCULATION \
		$(ENABLE_LOAD_CALCULATION)
	@printf '%-20s %s\n' TURBULENCE_WALL \
		$(ENABLE_TURBULENCE_WALL_TREATMENT)
	@printf '%-20s %s\n' MDF_PRESSURE_PROJ \
		$(ENABLE_MDF_PRESSURE_PROJECTION)
	@printf '%-20s %s\n' IBM_AWARE_PROJ \
		$(ENABLE_IBM_AWARE_PROJECTION)
	@printf '%-20s %s\n' ONE_WAY_DPM $(ENABLE_ONE_WAY_DPM)

check-openfoam:
	@test -n "$(WM_PROJECT_DIR)" || { \
		echo 'OpenFOAM is not loaded.' >&2; \
		exit 1; \
	}
	@command -v wmake >/dev/null 2>&1 || { \
		echo 'wmake is not available in PATH.' >&2; \
		exit 1; \
	}
	@command -v mpicxx >/dev/null 2>&1 || { \
		echo 'mpicxx is not available in PATH.' >&2; \
		exit 1; \
	}

help:
	@printf '%s\n' \
		'make              Build ichiBoundaryFoam.' \
		'make ibm          Clean and build the IBM solver.' \
		'make ibm EXE_NAME=mySolver  Build under another name.' \
		'make ibm EXE_PATH=/path/name  Select the full output file.' \
		'make ibm-aware EXE_NAME=ichiIbmAwareFoam  Build the experimental projection solver.' \
		'make ibm ENABLE_ONE_WAY_DPM=1  Add passive DPM tracking.' \
		'make ibm ENABLE_TURBULENCE_WALL_TREATMENT=0  Omit the immersed turbulence-wall model.' \
		'make pure-interfoam  Clean and build the pure baseline.' \
		'make marker-subiteration  Build the Constant--Pinelli prototype under a separate binary name.' \
		'make marker-operator-test  Build the no-PDE marker operator tester.' \
		'make clean-marker-operator-test  Clean only the tester objects.' \
		'make debug        Build with OpenFOAM debug flags.' \
		'make clean        Remove intermediate files.' \
		'make distclean    Also remove the executable.' \
		'make show-config  Print the build configuration.'

# vim: noexpandtab tabstop=4 shiftwidth=4
