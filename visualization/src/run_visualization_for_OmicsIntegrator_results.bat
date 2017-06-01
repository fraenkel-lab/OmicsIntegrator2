::
:: --------------------------------------------------------------
:: This script runs Visualization for the OmicsIntegrator results
:: --------------------------------------------------------------
::
SET FOREST_OUTPUT_DIR=%1
SET USER_SAMPLE=%2
SET USER_RESULT_DIR="../visualize_results_%USER_SAMPLE%"
mkdir %USER_RESULT_DIR%
::
::
:: Concatenate the graph_head with the graph_json.json (which is the output from Forest).
::
copy /B "graph_head.txt" + "%FOREST_OUTPUT_DIR%/graph_json.json" "%USER_RESULT_DIR%/graph.js"
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
py create_style_code_for_visualization.py %USER_RESULT_DIR% %FOREST_OUTPUT_DIR%
::
::
:: Make the HTML file that contains JavaScript code using the cytoscape.js library to be user-specific in that directory,
:: so the HTML can be linked to a domain name and the user can access the weblink for their result
::
copy visualize_OmicsIntegrator_results_with_CytoscapeJS_lib.htm %USER_RESULT_DIR%/visualize_OmicsIntegrator_results_with_CytoscapeJS_lib_%USER_SAMPLE%.htm
::
::
