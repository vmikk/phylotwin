#!/bin/bash

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT [-r RESOLUTION] [-t THREADS] [-m MEMORY] [-x TEMP_DIR]"
    echo "  -i INPUT          : Input directory containing Parquet files"
    echo "  -o OUTPUT         : Output Parquet file path"
    echo "  -r H3_RESOLUTION  : H3 resolution (e.g., 4)"
    echo "  -t THREADS        : Number of CPU threads to use (optional)"
    echo "  -m MEMORY         : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR       : Temporary directory path (optional)"
    echo "  -e EXT_DIR        : DuckDB extensions directory path (optional)"
    exit 1
}

INPUT=""
OUTPUT=""
H3_RESOLUTION=""

SQL_COMMAND+="

-- Install and load H3 extension
-- INSTALL h3 FROM community;
LOAD h3;

-- Create a table for WKT strings
CREATE TABLE polygons (
    wkt TEXT
);

-- Read WKT strings from input file
COPY polygons(wkt) FROM '${INPUT}' (
    DELIMITER '\t',
    HEADER FALSE
);

-- Convert WKT to H3 cells and save unique cells to CSV
COPY (
    SELECT DISTINCT
        UNNEST(h3_polygon_wkt_to_cells_string(wkt, ${H3_RESOLUTION})) AS h3_cell
    FROM polygons
    ORDER BY h3_cell
) TO '${OUTPUT}' (HEADER, DELIMITER ',');

-- -- For debugging: print number of cells per polygon
-- SELECT 
--     ROW_NUMBER() OVER () as polygon_id,
--     ARRAY_LENGTH(h3_polygon_wkt_to_cells_string(wkt, ${H3_RESOLUTION})) as cell_count
-- FROM polygons;

-- Clean up
DROP TABLE polygons;
"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"

