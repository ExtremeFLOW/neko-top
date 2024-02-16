# Define the type of driver to use.
# Options are:
#   - default: Clean neko driver
#   - user: Neko driver with user defined elements
#   - topopt: TopOpt driver
#   - custom: Custom driver
set(DRIVER_TYPE "default")

# If example require a custom driver, define the source file here, or store it
# as driver.f90 in the example directory.
# Otherwise leave commented out.
# set(DRIVER PATH/TO/CUSTOM/DRIVER.f90)

# Extra sources that need to be compiled into the driver.
# set(EXTRA_SOURCES "")

build_example(${DRIVER_TYPE})
