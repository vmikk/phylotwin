#!/bin/bash

SQL_COMMAND+="

-- Install and load H3 extension
-- INSTALL h3 FROM community;
LOAD h3;

-- Create a table for WKT strings
CREATE TABLE polygons (
    wkt TEXT
);

"

