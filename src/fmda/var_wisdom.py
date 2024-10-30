import numpy as np

_var_wisdom = {
    "dfm" : {
        "name" : "Dead Fuel moisture",
        "native_unit" : "gm",
        "colorbar" : "%",
        "colormap" : "jet_r",
        "norm_opt" : "boundary",
        "bounds" : [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3],
        "colors" : np.array([
            (156,22,27), (188,28,32), (217,45,43),
            (234,84,43), (245,137,56), (249,201,80),
            (215,225,95), (203,217,88), (114,190,75),
            (74,167,113), (60,150,120)
        ])/255.,
        "spacing" : "uniform",
        "scale" : [0.0, 0.3],
    },
    "lfm" : {
        "name" : "Live Fuel moisture",
        "native_unit" : "gm",
        "colorbar" : "%",
        "colormap" : "jet_r",
        "norm_opt" : "boundary",
        "bounds" : [0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.5, 2.0, 2.5, 3.0],
        "colors" : np.array([
            (156,22,27), (188,28,32), (217,45,43),
            (234,84,43), (245,137,56), (249,201,80),
            (215,225,95), (203,217,88), (114,190,75),
            (74,167,113), (60,150,120)
        ])/255.,
        "spacing" : "uniform",
        "scale" : [0.0, 3.0],
        "marker" : "^"
    },
    "EQUILd FM" : {
        "name" : "Drying equilibrium FM",
        "native_unit" : "-",
        "colorbar" : "-",
        "colormap" : "jet_r",
        "scale" : [0.0, 0.4]
    },
    "EQUILw FM" : {
        "name" : "Wetting equilibrium FM",
        "native_unit" : "-",
        "colorbar" : "-",
        "colormap" : "jet_r",
        "scale" : [0.0, 0.4]
    },
    "RH" : {
        "name" : "Relative humidity",
        "native_unit" : "%",
        "colorbar" : "%",
        "colormap" : "viridis_r",
        "scale" : [0.0, 100.0]
    },
    "T2" : {
        "name" : "Temperature at 2m",
        "native_unit" : "K",
        "colorbar" : "F",
        "colormap" : "jet",
        "scale" : [270.0, 320.0]
    },
    "PRECIP" : {
        "name" : "Precipitation",
        "native_unit" : "mm/h",
        "colorbar" : "mm/h",
        "colormap" : "jet_r",
        "transparent_values": [-np.inf, 0.2],
        "scale" : [0.0, 2.0]
    },
    "SNOWH" : {
        "name" : "Snow depth",
        "native_unit" : "m",
        "colorbar" : "mm",
        "colormap" : "tab20b_r",
        "transparent_values": [-np.inf, 0.002],
        "scale" : [0.0, 0.2]
    },
    "WINDSPD" : {
        "name" : "Wind Speed at 10 m",
        "native_unit" : "m/s",
        "colorbar" : "m/s",
        "colormap" : "jet",
        "bounds" : [0.0, 2.0, 5.0, 13.0, 18.0, 22.0, 27.0, 31.0],
        "colors" : np.array([
            (26,152,80), (102,189,99), (166,217,106),
            (217,239,139), (254,224,139), (253,174,97),
            (244,109,67), (215,48,39)
        ])/255.,
        "spacing" : "uniform",
        "scale" : [0.0, 32.0]
    },
    "WINDVEC" : {
        "name" : "Wind speed at 10 m",
        "components": ["U10", "V10"],
        "native_unit": "m/s",
        "scale": [0.0, 20.0]
    },
    "HGT" : {
        "name" : "Terrain height",
        "native_unit" : "m",
        "colorbar" : "m",
        "colormap" : "terrain",
        "scale" : [-86.0, 4500.0]
    },
    "SMOKE" : {
        "name" : "Smoke concentration",
        "native_unit" : "kg/m^3",
        "colorbar" : "ug/m^3",
        "colormap" : "rainbow",
        "norm_opt" : "boundary",
        "bounds" : [
            0, 1e-9, 2e-9, 4e-9, 6e-9, 8e-9, 
            1.2e-8, 1.6e-8, 2e-8, 2.5e-8, 3e-8, 
            4e-8, 6e-8, 1e-7, 2e-7
        ],
        "labels" : [
            "1", "2", "4", "6", "8", "12", "16",
            "20", "25", "30", "40", "60", "100",
            "200"
        ],
        "colors" : np.array([
            (255,255,255), (197,234,252), (148,210,240),
            (107,170,213), (72,149,176), (74,167,113),
            (114,190,75), (203,217,88), (249,201,80),
            (245,137,56), (234,84,43), (217,45,43),
            (188,28,32), (156,22,27), (147,32,205)
        ])/255.,
        "spacing" : "uniform",
        "transparent_values" : [-np.inf, 5e-10],
        "scale" : [0, 5e-7]
    },
    'FOSBERG' : {
        'name' : 'Fosberg Index',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'autumn_r',
        'scale' : [0, 45]
    },
    'HDW' : {
        'name' : 'Hot, Dry, & Windy Index',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'rainbow',
        'norm_opt' : 'boundary',
        'bounds' : [0,25,50,75,90,95,100],
        'colors' : np.array([
            (255,255,255), (254,239,180), (254, 200, 108),
            (243,147, 70), (198,104, 54), (135,  81,  56)
        ])/255.,
        'spacing' : 'uniform',
        'scale' : [0, 100]
    }
}

def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return _var_wisdom[var_name]
