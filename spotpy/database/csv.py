



from .base import database
import numpy as np
import io
import time
import sys
if sys.version_info[0] >= 3:
    str = str

class csv(database):
    """
    This class saves the process in the working storage. It can be used if
    safety matters.
    """

    def __init__(self, *args, **kwargs):
        # init base class
        super(csv, self).__init__(*args, **kwargs)
        # store init item only if dbinit
        if kwargs.get('dbappend', False) is False:
            print("* Database file '{}.csv' created.".format(self.dbname))
            # Create a open file, which needs to be closed after the sampling
            self.db = io.open(self.dbname + '.csv', 'w')
            # write header line
            self.db.write(str(','.join(self.header) + '\n'))
        else:
            print("* Appending to database file '{}.csv'.".format(self.dbname))
            # Continues writing file
            self.db = io.open(self.dbname + '.csv', 'a')

    def save(self, objectivefunction, parameterlist, simulations=None, chains=1):
        coll = (self.dim_dict['like'](objectivefunction) +
                self.dim_dict['par'](parameterlist) +
                self.dim_dict['simulation'](simulations) +
                [chains])
        # Apply rounding of floats
        coll = list(map(self.db_precision, coll))
        self.db.write(
            ','.join(map(str, coll)) + '\n')

        acttime = time.time()
        # Force writing to disc at least every two seconds
        if acttime - self.last_flush >= 2:
            self.db.flush()
            self.last_flush = time.time()

    def finalize(self):
        self.db.flush()  # Just to make sure that everything is written in file
        self.db.close()

    def getdata(self):
        data = np.genfromtxt(self.dbname + '.csv', delimiter=',', names=True)
        return data
