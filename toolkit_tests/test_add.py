import pytest
import copy
from pathlib import Path

import pandas as pd

import metobs_toolkit

libfolder = Path(str(Path(__file__).resolve())).parent.parent
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir

def get_demo_dataset():
    dataset = metobs_toolkit.Dataset()
    dataset.import_data_from_file(
        template_file=metobs_toolkit.demo_template,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        input_data_file=metobs_toolkit.demo_datafile,
    )
    return dataset

class TestAddMethods:
    def test_add_sensordata(self):
        # Get two SensorData objects with non-overlapping timestamps
        ds = get_demo_dataset()
        sta = ds.get_station("vlinder01")
        sd1 = copy.deepcopy(sta.get_sensor("temp"))

        #1. test without overlapping timestamps
        sd_other = copy.deepcopy(sta.get_sensor("temp"))
        sd_other.series.index = sd_other.series.index + pd.Timedelta('29D')
        sd_sum = sd1 + sd_other

        assert sd_sum.start_datetime == sd1.start_datetime
        assert sd_sum.end_datetime == sd_other.end_datetime
        #Note that there must be a gap (NAN's sandwiched, so do not check: len(sd_sum.df) == len(sd1.df) + len(sd2.df) )
        assert len(sd_sum.gaps) == 1
        n_records = int((sd_other.end_datetime - sd1.start_datetime)/sd1.freq)
        assert sd_sum.df.shape[0] == n_records+1

        #2. test with overlap in timestamps
        sd_other = copy.deepcopy(sta.get_sensor("temp"))
        sd_other.series.index = sd_other.series.index + pd.Timedelta('3D')
        sd_sum = sd1 + sd_other

        assert sd_sum.start_datetime == sd1.start_datetime
        assert sd_sum.end_datetime == sd_other.end_datetime
        #Note that there must be a gap (NAN's sandwiched, so do not check: len(sd_sum.df) == len(sd1.df) + len(sd2.df) )
        assert len(sd_sum.gaps) == 0
        n_records = int((sd_other.end_datetime - sd1.start_datetime)/sd1.freq)
        assert sd_sum.df.shape[0] == n_records+1

    def test_add_sensordata_different_obstypes(self):
        ds = get_demo_dataset()
        sta = ds.get_station("vlinder01")
        sd1 = copy.deepcopy(sta.get_sensor("temp"))
        sd_other = copy.deepcopy(sta.get_sensor("temp"))
        sd_other.series.index = sd_other.series.index + pd.Timedelta('3D')
        
        #1 test with incompatible obtype merge
        fake_obstype = metobs_toolkit.Obstype(obsname='temp',
                                              std_unit='degK', #makes for incompatibel merge! 
                                              description='blablabla')
        fake_obstype.original_unit='degK'
        fake_obstype.orginal_name="this_should be irrelevant"
        
        sd_other.obstype = fake_obstype
        # Test that adding SensorData with incompatible obstypes raises MetObsAdditionError
        from metobs_toolkit.backend_collection.errorclasses import MetObsAdditionError
        with pytest.raises(MetObsAdditionError):
            _ = sd1 + sd_other

    
        #2. test with compatible obstype merge
        
        fake_obstype = metobs_toolkit.Obstype(obsname='temp',
                                              std_unit='degC', #makes for compatibel merge! 
                                              description='blablabla')
        fake_obstype.original_unit='degK'
        fake_obstype.orginal_name="this_should be irrelevant"
        
        sd_other.obstype = fake_obstype
        
        sd_sum = sd1 + sd_other
        assert sd_sum.obstype.name == sd1.obstype.name
        assert sd_sum.obstype.std_unit == sd1.obstype.std_unit
        assert sd_sum.obstype.description != sd1.obstype.description #combination of both


    def test_add_sensordata_different_timezones(self):
        
        ds = get_demo_dataset()
        sta = ds.get_station("vlinder01")
        sd1 = copy.deepcopy(sta.get_sensor("temp"))
        sd_other = copy.deepcopy(sta.get_sensor("temp"))
        sd_other.series.index = sd_other.series.index + pd.Timedelta('3D')

        max_timestamp = sd_other.end_datetime

        #set timezone differntly
        # Make sd_other timestamps timezone-aware (assume original are naive/UTC), then convert to Asia/Shanghai
        sd_other.series.index = sd_other.series.index.tz_convert('Asia/Shanghai')
        sd_sum = sd1 + sd_other

        assert sd_sum.end_datetime == max_timestamp
        
        
    def test_add_sensordata_with_qc(self):

        ds = get_demo_dataset()
        sta = ds.get_station("vlinder01")
        sd_other = copy.deepcopy(sta.get_sensor("temp"))
        sd_other.series.index = sd_other.series.index + pd.Timedelta('30D')
       
        #apply qc before add
        sta.gross_value_check(target_obstype='temp', upper_threshold=24.2)
        sd1 = copy.deepcopy(sta.get_sensor("temp"))

        #summ
        sd_sum = sd1 + sd_other

        assert sd1.outliersdf.shape[0] !=0
        assert sd_sum.outliersdf.shape[0] == 0 #outliers must be removed


    def test_add_site_with_different_attributes(self):
        ds = get_demo_dataset()
        sta = ds.get_station("vlinder01")
        sta1 = copy.deepcopy(sta)
        sta1.get_landcover_fractions(buffers=[50])
        site1 = sta1.site

        sta2 = copy.deepcopy(sta)
        sta2.get_LCZ()
        sta2.get_landcover_fractions(buffers=[100,200])
        site2 = sta2.site
        
        # Add the two Site instances
        site_sum = site1 + site2
        # The merged site should keep the attributes from site1
        assert site_sum._stationname == site1._stationname
        assert site_sum.LCZ == site2.LCZ
        assert (len(site_sum._gee_buffered_fractions.keys()) == 3)


        #Incompatible merge test
        site1._lat=site1.lat + 0.0002
        from metobs_toolkit.backend_collection.errorclasses import MetObsAdditionError
        with pytest.raises(MetObsAdditionError):
            _ = site1 + site2
      

        

    def test_add_station_with_other_sensordata(self):
        ds = get_demo_dataset()
        sta_orig = copy.deepcopy(ds.get_station("vlinder01"))
        sta1 = copy.deepcopy(ds.get_station("vlinder01"))
        sta2 = copy.deepcopy(ds.get_station("vlinder01"))

        #Split sensordata
        sta1.obsdata = dict(list(sta1.obsdata.items())[:2])
        sta2.obsdata = dict(list(sta2.obsdata.items())[2:])
        sta_sum = sta1 + sta2
        assert_equality(to_check=sta_sum, solution=sta_orig)
        
        #split in time
        #Already tested on sensor level

        # add two different stations
        sta1 = copy.deepcopy(ds.get_station("vlinder01"))
        sta2 = copy.deepcopy(ds.get_station("vlinder02"))
        from metobs_toolkit.backend_collection.errorclasses import MetObsAdditionError
        with pytest.raises(MetObsAdditionError):
            _ = sta1 + sta2

        #conflict in coordinate
        sta1 = copy.deepcopy(ds.get_station("vlinder01"))
        sta2 = copy.deepcopy(ds.get_station("vlinder01"))
        sta2.site._lat = sta2.site.lat + 0.001
        from metobs_toolkit.backend_collection.errorclasses import MetObsAdditionError
        with pytest.raises(MetObsAdditionError):
            _ = sta1 + sta2


        
    def test_add_dataset(self):
        ds_orig = get_demo_dataset()
        ds1 = copy.deepcopy(ds_orig)
        ds2 = copy.deepcopy(ds_orig)

        dscomb = ds1+ds2
        for obs in dscomb.obstypes.values():
            trgobs = ds1.obstypes[obs.name]
            print(f'obs in comb: {obs} --> {type(obs)}')
            print(f'obs in ds1: {trgobs} --> {type(trgobs)}')
            # print(obs.get_info())
            print(trgobs.get_info())
            print(obs == trgobs)

        # assert_equality(to_check=ds1+ds2, solution=ds_orig)

        #1 Split in two sets by station
        # ds1.stations = ds1.stations[:16] 
        # ds2.stations = ds2.stations[16:]
        # assert_equality(to_check=ds1+ds2, solution=ds_orig)


        #2 merge metadata only (with metata) and with data





if __name__ == "__main__":
    t = TestAddMethods()
    # t.test_add_sensordata()
    # t.test_add_station()
    # t.test_add_site()
    # t.test_add_modeltimeseries()
    t.test_add_dataset()
