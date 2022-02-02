Setup:

In an mrb area, clone this repository into the srcs/ area. Then run `mrb
uc` and build the mrb area is normally

Generate events (currently up to detsim) -- run in an area where you
have already run `mrbslp`:

`lar -c prodsingle_larbox.fcl -n 10`

Display events:

`lar -c evd_larbox.fcl -s /path/to/file.root`

