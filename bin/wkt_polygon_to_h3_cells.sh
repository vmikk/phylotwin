#!/bin/bash

## Convert WKT polygons to H3 cells

## NB. 
## - Input is a text file with polygons specified in WKT format
##    Each line is a separate polygon
##    (MULTIPOLYGONs and POINTs are currently not supported)
## - Expected CRS is EPSG:4326

## Notes on MULTIPOLYGON:
# - workaround ST_Dump(geom).unnest(recursive := true) + flaten(list())...
#   https://github.com/isaacbrodsky/h3-duckdb/issues/118



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

## Initialize variables
INPUT=""
OUTPUT=""
H3_RESOLUTION=""
THREADS=""
MEMORY=""
TEMP_DIR=""
EXT_DIR=""

## Parse command-line options
while getopts "i:o:r:t:m:x:e" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) H3_RESOLUTION="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        e) EXT_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

echo "Validating input parameters"

## Validate input parameters
if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$H3_RESOLUTION" ]]; then
    echo -e "Error: Missing required parameters (INPUT, OUTPUT, H3_RESOLUTION)!\n"
    usage
fi

## H3 resolution should be an integer between 0 and 15
if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

echo "Counting user-supplied WKT polygons"

## Count number of polygons
POLYGON_COUNT=$(grep -c "^POLYGON" "${INPUT}")
NUM_LINES=$(wc -l < "${INPUT}")
echo "..Number of lines in input file: $NUM_LINES"
echo "..Number of polygons detected: $POLYGON_COUNT"
if [ "$POLYGON_COUNT" -ne "$NUM_LINES" ]; then
    echo -e "WARNING: Number of polygons detected does not match number of lines in input file!\n"
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "..Input:  $INPUT"
echo "..Output: $OUTPUT"
echo "..H3 resolution: $H3_RESOLUTION"

if [[ -n "$THREADS" ]]; then
    echo "..Threads: $THREADS"
fi
if [[ -n "$MEMORY" ]]; then
    echo "..Memory: $MEMORY"
fi
if [[ -n "$TEMP_DIR" ]]; then
    echo "..Temp directory: $TEMP_DIR"
fi
if [[ -n "$EXT_DIR" ]]; then
    echo "..DuckDB extensions directory: $EXT_DIR"
fi

## Start the SQL command
echo -e "\nPreparing SQL command"

SQL_COMMAND=""

## Add configuration settings (if provided)
if [[ -n "$THREADS" ]]; then
    SQL_COMMAND+="
SET threads TO ${THREADS};
"
fi

if [[ -n "$MEMORY" ]]; then
    SQL_COMMAND+="
SET memory_limit = '${MEMORY}';
"
fi

if [[ -n "$TEMP_DIR" ]]; then
    SQL_COMMAND+="
PRAGMA temp_directory='${TEMP_DIR}';
"
fi

if [[ -n "$EXT_DIR" ]]; then
    SQL_COMMAND+="
SET extension_directory='${EXT_DIR}';
"
fi

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

