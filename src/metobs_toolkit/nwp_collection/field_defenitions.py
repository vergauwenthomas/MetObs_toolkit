
from metobs_toolkit.obstypes import tlk_obstypes, Obstype, ModelObstype

default_SFX_fields = [
    
    ModelObstype(obstype=tlk_obstypes['temp'],
                                 model_unit='degK',
                                 model_band="SFX.T2M_ISBA"),
    ModelObstype(obstype=tlk_obstypes['temp'],
                                 model_unit='degK',
                                 model_band="SFX.T2M_TEB"),
    ModelObstype(obstype=tlk_obstypes['temp'],
                                 model_unit='degK',
                                 model_band="SFX.T2M"),

    ModelObstype(obstype=tlk_obstypes['wind_speed'],
                                 model_unit='m/s',
                                 model_band="SFX.ZON10M"),
]


default_AROME_fields = [
    ModelObstype(obstype=tlk_obstypes['temp'],
                model_unit='degK',
                model_band="CLSTEMPERATURE"),
]
