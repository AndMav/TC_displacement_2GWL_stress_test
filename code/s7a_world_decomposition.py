"""
Merges the decomposition at every exposure point of all countries into one table for the world maps: the dominant
driver and its value at every point with a dominant driver, at p10, median and p90.

Input:  {country}_EAD_exp_point_decomposition.csv         (decomposition at every exposure point, every country)
Output: world_points.csv                                  (country, latitude, longitude, stat, dominant_driver,
                                                           dominant_value)
        (folders defined in config.py)
"""

import pandas as pd
from config import DECOMPOSITION_DIR, WORLD_DATA_DIR, COUNTRIES

'''
 ----------------------------------------------------------------------------------------
     I. Constants
 ----------------------------------------------------------------------------------------
'''

# Columns of the world table
COLUMNS = ['country', 'latitude', 'longitude', 'stat', 'dominant_driver', 'dominant_value']


'''
 ----------------------------------------------------------------------------------------
     II. World table
 ----------------------------------------------------------------------------------------
'''

def world_points():
    """Exposure points of all countries, for every statistic, with their dominant driver.
    Output: one dataframe, one row per exposure point and statistic (all countries stacked), with the columns
        country, latitude, longitude, stat, dominant_driver, dominant_value
    """

    tables = []

    for country in COUNTRIES:

        path = DECOMPOSITION_DIR / country / f"{country}_EAD_exp_point_decomposition.csv"
        df = pd.read_csv(path)
        df = df[df['dominant_driver'] != 'NO_DATA']    # points without a dominant driver are not drawn
        table = df[COLUMNS]
        tables.append(table)

    return pd.concat(tables, ignore_index=True)

'''
 ----------------------------------------------------------------------------------------
     III. Computation and saving of the table
 ----------------------------------------------------------------------------------------
'''

def main():
    
    WORLD_DATA_DIR.mkdir(parents=True, exist_ok=True)
    points = world_points()
    points.to_csv(WORLD_DATA_DIR / "world_points.csv", index=False)


    print(f"Done -> {WORLD_DATA_DIR / 'world_points.csv'}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
