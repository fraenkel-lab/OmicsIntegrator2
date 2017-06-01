::
:: --------------------------------------------------------------
:: This script runs Visualization for the user saved JSON file
:: --------------------------------------------------------------
::
SET USER_SAVED_JSON_DIR=%1
SET USER_SAMPLE=%2
SET USER_RESULT_DIR="../visualize_saved_JSON_%USER_SAMPLE%"
mkdir %USER_RESULT_DIR%
::
::
:: Concatenate the graph_head with the graph_json.json (which is the user saved JSON file).
::
copy /B "graph_head.txt" + "%USER_SAVED_JSON_DIR%/graph_json.json" "%USER_RESULT_DIR%/graph.js"
::
::
:: Copy the style_master to the visualize_results_rundir/
::
copy "style_master.js" "%USER_RESULT_DIR%/style.js"
::
::
:: Run Python code to extract min and max values from graph_json.json, to create the style.js code, 
:: and to process the graph.js.
::
py create_style_code_for_visualization.py %USER_RESULT_DIR% %USER_SAVED_JSON_DIR%
::
::
:: Make the HTML file that contains JavaScript code using the cytoscape.js library to be user-specific in that directory,
:: so the HTML can be linked to a domain name and the user can access the weblink for their result
::
copy visualize_user_saved_JSON_file_with_CytoscapeJS_lib.htm %USER_RESULT_DIR%/visualize_user_saved_JSON_file_with_CytoscapeJS_lib_%USER_SAMPLE%.htm
::
::
