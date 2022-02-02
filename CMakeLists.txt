# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


cmake_minimum_required (VERSION 3.12...3.18 FATAL_ERROR)

project(larbox LANGUAGES CXX)

message(STATUS "\n\n  ==========================   ${PROJECT_NAME}   ==========================")

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules $ENV{ART_DIR}/Modules )

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-undefined-var-template
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product(art)
find_ups_product(art_root_io)
find_ups_product(gsl)
find_ups_product(nusimdata)
find_ups_product(nug4)
find_ups_product(lardataobj)
find_ups_product(larcorealg)
find_ups_product(larcore)
find_ups_product(larreco)
find_ups_product(larpandora)
find_ups_product(larsoft v09_00_00)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(larbox)
add_subdirectory(ups)
add_subdirectory(fcl)
add_subdirectory(data)

# packaging utility
include(UseCPack)
