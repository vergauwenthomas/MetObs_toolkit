# -*- coding: utf-8 -*-

""" Test the dependancies in virtual environments.
 Just run 'nox' in the terminal to run it locally. 

(note that only installed python versions are runned!)
 """
import nox


@nox.session(python=["3.7", "3.8", "3.9", "3.10"])
def test(session):
    session.install('pandas == 1.3.0',
                    'numpy == 1.17.3',
                    'matplotlib == 3.0.0',
                    'mysql-connector-python >= 8.0.6',
                    'geopandas == 0.9.0',
                    'rasterstats == 0.14.0',
                    'mapclassify == 2.4.0',
                    'openpyxl == 3.0.0',
                    'xlrd >= 1.0.0')
    
    
    
    session.run("python", "tests/push_test/basic_import_test.py")
    session.run("python", "tests/push_test/IO_test.py")
    session.run("python", "tests/push_test/plot_test.py")
    session.run("python", "tests/push_test/qc_test.py")
