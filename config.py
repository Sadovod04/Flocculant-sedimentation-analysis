"""Shared thickener / feedwell constants.

These used to be read from a binary ``database.db`` SQLite file that is not
reproducible (and had become corrupted). Keeping them here as plain Python makes
every script self-contained and the numbers reviewable in a diff.

Units follow the source coursework: flows in m3/h, densities in kg/m3,
viscosity in Pa*s, lengths in m.
"""

# --- thickener geometry ---
CONE_HEIGHT = 1.35        # height of the conical bottom, m
CYLINDER_HEIGHT = 1.0     # height of the cylindrical section below the feedwell, m
DIAMETER = 30.0           # thickener diameter, m

# --- volumetric flows ---
FEED_FLOW = 350.0         # underflow feed slurry, m3/h  (Qufeed)
DILUTION_FLOW = 20.0      # feedwell dilution, m3/h       (Qinj)
UNDERFLOW = 90.0          # underflow withdrawal, m3/h    (Qunderfl)

# --- phase properties ---
FEED_SOLID_FRACTION = 0.0159   # volumetric solids fraction in the feed (Fifeed)
SOLID_DENSITY = 3200.0         # particle density, kg/m3  (psolid)
LIQUOR_DENSITY = 1240.0        # liquor density, kg/m3    (pfluid)
WATER_DENSITY = 1020.0         # dilution water density, kg/m3
LIQUOR_VISCOSITY = 0.0021      # liquor dynamic viscosity, Pa*s (muliqour)

# --- flocculant ---
FLOCCULANT_FLOW_PCT = 2.0      # flocculant addition, %   (Qfloc)
FLOCCULANT_CONC = 0.005        # flocculant working solution concentration (Cfloc_w)

GRAVITY = 9.81

# File the population-balance model writes its result to, and the thickener
# profile model reads the mean floc diameter from.
MEAN_DIAMETER_FILE = "data/mean_floc_diameter.json"
