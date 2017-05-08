#!/bin/bash
# --------------------------------------------------------------
# This script runs Visualization for the user saved JSON file
# --------------------------------------------------------------
#
USER_SAVED_JSON_DIR=$1
USER_SAMPLE=$2
USER_RESULT_DIR="../visualize_saved_JSON_$USER_SAMPLE"
#
# Make a visualize_results_rundir/ after clearing out the existing one (if exists). 
# This run dir is for creating the final graph.js and style.js files of a specific network, 
# which will be used to visualize the network using the cytoscape.js library.
#
rm -r ../visualize_results_rundir
mkdir ../visualize_results_rundir
#
#
# Concatenate the graph_head with the graph_json.json (which is the user saved JSON file).
#
cp $USER_SAVED_JSON_DIR/graph_json.json .
cat graph_head.txt graph_json.json > ../visualize_results_rundir/graph.js
#
#
# Copy the style_master to the visualize_results_rundir/
#
cp style_master.js ../visualize_results_rundir/style.js
#
#
# Run Python code to extract min and max values from graph_json.json, to create the style.js code, 
# and to process the graph.js.
#
python create_style_code_for_visualization.py
rm graph_json.json
#
#
# Now copy graph.js and style.js to a user-specific directory
# and make the HTML file that contains JavaScript code using the cytoscape.js library to be user-specific in that directory,
# so the HTML can be linked to a domain name and the user can access the weblink for their result
#
mkdir $USER_RESULT_DIR
cp ../visualize_results_rundir/graph.js $USER_RESULT_DIR
cp ../visualize_results_rundir/style.js $USER_RESULT_DIR
cp visualize_user_saved_JSON_file_with_CytoscapeJS_lib.htm $USER_RESULT_DIR/visualize_user_saved_JSON_file_with_CytoscapeJS_lib_$USER_SAMPLE.htm
#
